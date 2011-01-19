//
// Copyright (C) 2010 Helene Lageber
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, see <http://www.gnu.org/licenses/>.
//

#include "BGPRouting.h"
#include "RoutingTableAccess.h"
#include "OSPFRouting.h"
#include "BGPSession.h"

Define_Module(BGPRouting);

bool            OSPFExist(IRoutingTable* rtTable);
simtime_t*      loadTimerConfig(cXMLElementList& timerConfig);
BGP::ASID       findMyAS(cXMLElementList& ASList, IRoutingTable* rtTable, unsigned char* routerPosition);
unsigned char   ASLoopDetection(BGP::RoutingTableEntry* entry, BGP::ASID myAS);
BGP::SessionID  findIdFromPeerAddr(std::map<BGP::SessionID, BGPSession*> sessions, IPAddress peerAddr);
int             isInIPTable(IRoutingTable* rtTable, IPAddress addr);
BGP::SessionID  findIdFromSocketConnId(std::map<BGP::SessionID, BGPSession*> sessions, int connId);
unsigned int    calculateStartDelay(int rtListSize, unsigned char rtPosition, unsigned char rtPeerPosition);

BGPRouting::~BGPRouting(void)
{
    for (std::map<BGP::SessionID, BGPSession*>::iterator sessionIterator = _BGPSessions.begin();
        sessionIterator != _BGPSessions.end(); sessionIterator ++)
    {
        (*sessionIterator).second->~BGPSession();
    }
    _BGPRoutingTable.erase(_BGPRoutingTable.begin(), _BGPRoutingTable.end());
    _prefixListIN.erase(_prefixListIN.begin(), _prefixListIN.end());
    _prefixListOUT.erase(_prefixListOUT.begin(), _prefixListOUT.end());
}

void BGPRouting::initialize(int stage)
{
    if (stage==4)
    {
        _rt = RoutingTableAccess().get();
        _inft = InterfaceTableAccess().get();

        // read BGP configuration
        const char *fileName = par ("bgpConfigFile");
        if (*fileName == 0 || !loadConfigFromXML (fileName))
        {
            error ("Error reading BGP configuration from file %s", fileName);
        }

        createWatch("myAutonomousSystem", _myAS);
        WATCH_PTRVECTOR(_BGPRoutingTable);
    }
}

void BGPRouting::handleMessage(cMessage *msg)
{
    if (msg->isSelfMessage()) //BGP level
    {
        handleTimer(msg);
    }
    else if (!strcmp(msg->getArrivalGate()->getName(), "fromTCP")) //TCP level
    {
        processMessageFromTCP(msg);
    }
    else
    {
        delete msg;
    }
}

void BGPRouting::handleTimer(cMessage *timer)
{
    BGPSession* pSession = (BGPSession*)timer->getContextPointer();
    if (pSession)
    {
        switch (timer->getKind())
        {
            case BGP::START_EVENT_KIND:
                EV << "Processing Start Event" << std::endl;
                pSession->getFSM()->ManualStart();
                break;
            case BGP::CONNECT_RETRY_KIND:
                EV << "Expiring Connect Retry Timer" << std::endl;
                pSession->getFSM()->ConnectRetryTimer_Expires();
                break;
            case BGP::HOLD_TIME_KIND:
                EV << "Expiring Hold Timer" << std::endl;
                pSession->getFSM()->HoldTimer_Expires();
                break;
            case BGP::KEEP_ALIVE_KIND:
                EV << "Expiring Keep Alive timer" << std::endl;
                pSession->getFSM()->KeepaliveTimer_Expires();
                break;
            default :
                error("Invalid timer kind %d",timer->getKind());
        }
    }
}

void BGPRouting::finish()
{
    unsigned int statTab[BGP::NB_STATS] = {0,0,0,0,0,0};
    for (std::map<BGP::SessionID, BGPSession*>::iterator sessionIterator = _BGPSessions.begin(); sessionIterator != _BGPSessions.end(); sessionIterator ++)
    {
        (*sessionIterator).second->getStatistics(statTab);
    }
    recordScalar("OPENMsgSent", statTab[0]);
    recordScalar("OPENMsgRecv", statTab[1]);
    recordScalar("KeepAliveMsgSent", statTab[2]);
    recordScalar("KeepAliveMsgRcv",statTab[3]);
    recordScalar("UpdateMsgSent", statTab[4]);
    recordScalar("UpdateMsgRcv", statTab[5]);
}

void BGPRouting::listenConnectionFromPeer(BGP::SessionID sessionID)
{
    if (_BGPSessions[sessionID]->getSocketListen()->getState() == TCPSocket::CLOSED)
    {
        //session StartDelayTime error, it's anormal that listenSocket is closed.
        _socketMap.removeSocket(_BGPSessions[sessionID]->getSocketListen());
        _BGPSessions[sessionID]->getSocketListen()->abort();
        _BGPSessions[sessionID]->getSocketListen()->renewSocket();
    }
    if (_BGPSessions[sessionID]->getSocketListen()->getState() != TCPSocket::LISTENING)
    {
        _BGPSessions[sessionID]->getSocketListen()->setOutputGate(gate("toTCP"));
        _BGPSessions[sessionID]->getSocketListen()->readDataTransferModePar(*this);
        _BGPSessions[sessionID]->getSocketListen()->bind(BGP::TCP_PORT);
        _BGPSessions[sessionID]->getSocketListen()->listen();
    }
}

void BGPRouting::openTCPConnectionToPeer(BGP::SessionID sessionID)
{
    InterfaceEntry* intfEntry = _BGPSessions[sessionID]->getLinkIntf();
    TCPSocket *socket = _BGPSessions[sessionID]->getSocket();
    if (socket->getState() != TCPSocket::NOT_BOUND)
    {
        _socketMap.removeSocket(socket);
        socket->abort();
        socket->renewSocket();
    }
    socket->setCallbackObject(this, (void*)sessionID);
    socket->setOutputGate(gate("toTCP"));
    socket->readDataTransferModePar(*this);
    socket->bind(intfEntry->ipv4Data()->getIPAddress(),0);
    _socketMap.addSocket(socket);

    socket->connect(_BGPSessions[sessionID]->getPeerAddr(), BGP::TCP_PORT);
}

void BGPRouting::processMessageFromTCP(cMessage *msg)
{
    TCPSocket *socket = _socketMap.findSocketFor(msg);
    if (!socket)
    {
        socket = new TCPSocket(msg);
        socket->readDataTransferModePar(*this);
        socket->setOutputGate(gate("toTCP"));
        IPAddress peerAddr = socket->getRemoteAddress().get4();
        BGP::SessionID i = findIdFromPeerAddr(_BGPSessions, peerAddr);
        if (i==-1)
        {
            socket->close();
            delete socket;
            delete msg;
            return;
        }
        socket->setCallbackObject(this, (void *)i);

        _socketMap.addSocket(socket);
        _BGPSessions[i]->getSocket()->abort();
        _BGPSessions[i]->setSocket(socket);
    }

    socket->processMessage(msg);
}

void BGPRouting::socketEstablished(int connId, void *yourPtr)
{
    _currSessionId = findIdFromSocketConnId(_BGPSessions, connId);
    if (_currSessionId == -1)
    {
        error("socket id=%d is not established",connId);
    }

    //if it's an IGP Session, TCPConnectionConfirmed only if all EGP Sessions established
    if (_BGPSessions[_currSessionId]->getType() == BGP::IGP &&
        this->findNextSession(BGP::EGP) != -1)
    {
        _BGPSessions[_currSessionId]->getFSM()->TcpConnectionFails();
    }
    else
    {
        _BGPSessions[_currSessionId]->getFSM()->TcpConnectionConfirmed();
        _BGPSessions[_currSessionId]->getSocketListen()->abort();
    }

    if (_BGPSessions[_currSessionId]->getSocketListen()->getConnectionId() != connId &&
        _BGPSessions[_currSessionId]->getType() == BGP::EGP &&
        this->findNextSession(BGP::EGP) != -1 )
    {
        _BGPSessions[_currSessionId]->getSocketListen()->abort();
    }
}

void BGPRouting::socketDataArrived(int connId, void *yourPtr, cPacket *msg, bool urgent)
{
    _currSessionId = findIdFromSocketConnId(_BGPSessions, connId);
    if (_currSessionId != -1)
    {
        BGPHeader* ptrHdr = check_and_cast<BGPHeader*>(msg);
        switch(ptrHdr->getType())
        {
            case BGP_OPEN:
                //BGPOpenMessage* ptrMsg = check_and_cast<BGPOpenMessage*>(msg);
                processMessage(*check_and_cast<BGPOpenMessage*>(msg));
                break;
            case BGP_KEEPALIVE:
                processMessage(*check_and_cast<BGPKeepAliveMessage*>(msg));
                break;
            case BGP_UPDATE:
                processMessage(*check_and_cast<BGPUpdateMessage*>(msg));
                break;
            default:
                error("Invalid BGP message type %d",ptrHdr->getType());
        }
    }
    delete msg;
}

void BGPRouting::socketFailure(int connId, void *yourPtr, int code)
{
    _currSessionId = findIdFromSocketConnId(_BGPSessions, connId);
    if (_currSessionId != -1)
    {
        _BGPSessions[_currSessionId]->getFSM()->TcpConnectionFails();
    }
}

void BGPRouting::processMessage(const BGPOpenMessage& msg)
{
    EV << "Processing BGP OPEN message" << std::endl;
    _BGPSessions[_currSessionId]->getFSM()->OpenMsgEvent();
}

void BGPRouting::processMessage(const BGPKeepAliveMessage& msg)
{
    EV << "Processing BGP Keep Alive message" << std::endl;
    _BGPSessions[_currSessionId]->getFSM()->KeepAliveMsgEvent();
}

void BGPRouting::processMessage(const BGPUpdateMessage& msg)
{
    EV << "Processing BGP Update message" << std::endl;
    _BGPSessions[_currSessionId]->getFSM()->UpdateMsgEvent();

    unsigned char               decisionProcessResult;
    IPAddress                   netMask(IPAddress::ALLONES_ADDRESS);
    BGP::RoutingTableEntry*     entry           = new BGP::RoutingTableEntry();
    const unsigned char         length          = msg.getNLRI().length;
    unsigned int                ASValueCount    = msg.getPathAttributeList(0).getAsPath(0).getValue(0).getAsValueArraySize();

    entry->setDestinationID(msg.getNLRI().prefix);
    netMask.keepFirstBits(32-length);
    entry->setAddressMask(netMask);
    for (unsigned int j=0; j< ASValueCount; j++)
    {
        entry->addAS(msg.getPathAttributeList(0).getAsPath(0).getValue(0).getAsValue(j));
    }

    decisionProcessResult = ASLoopDetection(entry, _myAS);

    if (decisionProcessResult == BGP::ASLOOP_NO_DETECTED)
    {
        // RFC 4271, 9.1.  Decision Process
        decisionProcessResult = decisionProcess(msg, entry, _currSessionId);
        //RFC 4271, 9.2.  Update-Send Process
        if (decisionProcessResult != 0)
        {
            updateSendProcess(decisionProcessResult, _currSessionId, entry);
        }
    }
}

unsigned char BGPRouting::decisionProcess(const BGPUpdateMessage& msg, BGP::RoutingTableEntry* entry, BGP::SessionID sessionIndex)
{
    //Don't add the route if it exists in PrefixListINTable or in ASListINTable
    if (isInTable(_prefixListIN, entry) != -1 || isInASList(_ASListIN, entry))
    {
        return 0;
    }

    /*If the AS_PATH attribute of a BGP route contains an AS loop, the BGP
    route should be excluded from the decision process. */
    entry->setPathType(msg.getPathAttributeList(0).getOrigin().getValue());
    entry->setNextHop(msg.getPathAttributeList(0).getNextHop().getValue());

    //if the route already exist in BGP routing table, tieBreakingProcess();
    //(RFC 4271: 9.1.2.2 Breaking Ties)
    unsigned long BGPindex = isInTable(_BGPRoutingTable, entry);
    if (BGPindex != -1)
    {
        if (tieBreakingProcess(_BGPRoutingTable[BGPindex], entry))
        {
            return 0;
        }
        else
        {
            entry->setInterface(_BGPSessions[sessionIndex]->getLinkIntf());
            _BGPRoutingTable.push_back(entry);
            _rt->addRoute(entry);
            return BGP::ROUTE_DESTINATION_CHANGED;
        }
    }

    //Don't add the route if it exists in IP routing table except if the msg come from IGP session
    int indexIP = isInIPTable(_rt, entry->getHost());
    if (indexIP != -1 && _rt->getRoute(indexIP)->getSource() != IPRoute::BGP )
    {
        if (_BGPSessions[sessionIndex]->getType() != BGP::IGP )
        {
            return 0;
        }
        else
        {
            IPRoute* newEntry = new IPRoute;
            newEntry->setHost(_rt->getRoute(indexIP)->getHost());
            newEntry->setNetmask(_rt->getRoute(indexIP)->getNetmask());
            newEntry->setGateway(_rt->getRoute(indexIP)->getGateway());
            newEntry->setInterface(_rt->getRoute(indexIP)->getInterface());
            newEntry->setSource(IPRoute::BGP);
            newEntry->setType(_rt->getRoute(indexIP)->getType());
            _rt->deleteRoute(_rt->getRoute(indexIP));
            _rt->addRoute(newEntry);
        }
    }

    entry->setInterface(_BGPSessions[sessionIndex]->getLinkIntf());
    _BGPRoutingTable.push_back(entry);

    if (_BGPSessions[sessionIndex]->getType() == BGP::EGP)
    {
        std::string entryh = entry->getHost().str();
        std::string entryn = entry->getNetmask().str();
        _rt->addRoute(entry);
        //insertExternalRoute on OSPF ExternalRoutingTable if OSPF exist on this BGP router
        if (OSPFExist(_rt))
        {
            OSPF::IPv4AddressRange  OSPFnetAddr;
            OSPFnetAddr.address = ipv4AddressFromULong(entry->getHost().getInt());
            OSPFnetAddr.mask = ipv4AddressFromULong(entry->getNetmask().getInt());
            OSPFRouting* ospf = OSPFRoutingAccess().getIfExists();
            ospf->insertExternalRoute(entry->getInterfaceName(), OSPFnetAddr);
            simulation.setContext(this);
        }
    }
    return BGP::NEW_ROUTE_ADDED;
}

bool BGPRouting::tieBreakingProcess(BGP::RoutingTableEntry* oldEntry, BGP::RoutingTableEntry* entry)
{
    /*a) Remove from consideration all routes that are not tied for
         having the smallest number of AS numbers present in their
         AS_PATH attributes.*/
    if (entry->getASCount() < oldEntry->getASCount())
    {
        deleteBGPRoutingEntry(oldEntry);
        return false;
    }

    /* b) Remove from consideration all routes that are not tied for
         having the lowest Origin number in their Origin attribute.*/
    if (entry->getPathType() < oldEntry->getPathType())
    {
        deleteBGPRoutingEntry(oldEntry);
        return false;
    }
    return true;
}

void BGPRouting::updateSendProcess(const unsigned char type, BGP::SessionID sessionIndex, BGP::RoutingTableEntry* entry)
{
    //Don't send the update Message if the route exists in listOUTTable
    //SESSION = EGP : send an update message to all BGP Peer (EGP && IGP)
    //if it is not the currentSession and if the session is already established
    //SESSION = IGP : send an update message to External BGP Peer (EGP) only
    //if it is not the currentSession and if the session is already established
    for (std::map<BGP::SessionID, BGPSession*>::iterator sessionIt = _BGPSessions.begin();
        sessionIt != _BGPSessions.end(); sessionIt ++)
    {
        if (isInTable(_prefixListOUT, entry) != -1  || isInASList(_ASListOUT, entry)      ||
            ((*sessionIt).first == sessionIndex && type != BGP::NEW_SESSION_ESTABLISHED ) ||
            (type == BGP::NEW_SESSION_ESTABLISHED && (*sessionIt).first != sessionIndex ) ||
            !(*sessionIt).second->isEstablished() )
        {
            continue;
        }
        if ((_BGPSessions[sessionIndex]->getType()==BGP::IGP && (*sessionIt).second->getType()==BGP::EGP ) ||
            _BGPSessions[sessionIndex]->getType() == BGP::EGP   ||
            type == BGP::ROUTE_DESTINATION_CHANGED              ||
            type == BGP::NEW_SESSION_ESTABLISHED )
        {
            BGPUpdateNLRI                   NLRI;
            BGPUpdatePathAttributeList  content;

            unsigned int nbAS = entry->getASCount();
            content.setAsPathArraySize(1);
            content.getAsPath(0).setValueArraySize(1);
            content.getAsPath(0).getValue(0).setType(BGP::AS_SEQUENCE);
            //RFC 4271 : set My AS in first position if it is not already
            if (entry->getAS(0) != _myAS)
            {
                content.getAsPath(0).getValue(0).setAsValueArraySize(nbAS+1);
                content.getAsPath(0).getValue(0).setLength(1);
                content.getAsPath(0).getValue(0).setAsValue(0, _myAS);
                for (unsigned int j = 1; j < nbAS+1; j++)
                {
                    content.getAsPath(0).getValue(0).setAsValue(j, entry->getAS(j-1));
                }
            }
            else
            {
                content.getAsPath(0).getValue(0).setAsValueArraySize(nbAS);
                content.getAsPath(0).getValue(0).setLength(1);
                for (unsigned int j = 0; j < nbAS; j++)
                {
                    content.getAsPath(0).getValue(0).setAsValue(j, entry->getAS(j));
                }
            }

            InterfaceEntry*  iftEntry = (*sessionIt).second->getLinkIntf();
            content.getOrigin().setValue((*sessionIt).second->getType());
            content.getNextHop().setValue(iftEntry->ipv4Data()->getIPAddress());
            IPAddress netMask = entry->getNetmask();
            NLRI.prefix = entry->getHost().doAnd(netMask);
            NLRI.length = (unsigned char) netMask.getNetmaskLength();
            {
                BGPUpdateMessage* updateMsg = new BGPUpdateMessage("BGPUpdate");
                updateMsg->setPathAttributeListArraySize(1);
                updateMsg->setPathAttributeList(content);
                updateMsg->setNLRI(NLRI);
                (*sessionIt).second->getSocket()->send(updateMsg);
                (*sessionIt).second->addUpdateMsgSent();
            }
        }
    }
}

bool BGPRouting::checkExternalRoute(const IPRoute* route)
{
    OSPF::IPv4Address OSPFRoute;
    OSPFRoute = ipv4AddressFromULong(route->getHost().getInt());
    OSPFRouting* ospf = OSPFRoutingAccess().getIfExists();
    bool returnValue = ospf->checkExternalRoute(OSPFRoute);
    simulation.setContext(this);
    return returnValue;
}

void loadTimerConfig(cXMLElementList& timerConfig, simtime_t* delayTab)
{
    for (cXMLElementList::iterator timerElemIt = timerConfig.begin(); timerElemIt != timerConfig.end(); timerElemIt++)
    {
        std::string nodeName = (*timerElemIt)->getTagName();
        if (nodeName == "connectRetryTime")
        {
            delayTab[0] = (double)atoi((*timerElemIt)->getNodeValue());
        }
        if (nodeName == "holdTime")
        {
            delayTab[1] = (double)atoi((*timerElemIt)->getNodeValue());
        }
        if (nodeName == "keepAliveTime")
        {
            delayTab[2] = (double)atoi((*timerElemIt)->getNodeValue());
        }
        if (nodeName == "startDelay")
        {
            delayTab[3] = (double)atoi((*timerElemIt)->getNodeValue());
        }
    }
}

BGP::ASID findMyAS(cXMLElementList& ASList, IRoutingTable* rtTable, unsigned char* routerPositionPtr)
{
    const char* routerNode;
    BGP::ASID myAS = 0;
    for (cXMLElementList::iterator ASListIt = ASList.begin(); ASListIt != ASList.end() && myAS == 0; ASListIt++)
    {
        cXMLElementList routerList = (*ASListIt)->getChildrenByTagName("Router");
        *routerPositionPtr = 1;
        for (cXMLElementList::iterator routerListIt = routerList.begin(); routerListIt != routerList.end(); routerListIt++)
        {
            routerNode = (*routerListIt)->getAttribute("interAddr");
            IPAddress routerAddr = IPAddress(routerNode);
            if (isInIPTable(rtTable, routerAddr) == -1)
            {
                *routerPositionPtr +=1;
                continue;
            }
            myAS = atoi((*routerListIt)->getParentNode()->getAttribute("id"));
            return myAS;
        }
    }

    return 0;
}

void BGPRouting::loadSessionConfig(cXMLElementList& sessionList, simtime_t* delayTab)
{
    simtime_t saveStartDelay = delayTab[3];
    for (cXMLElementList::iterator sessionListIt = sessionList.begin(); sessionListIt != sessionList.end(); sessionListIt++, delayTab[3] = saveStartDelay)
    {
        const char* exterAddr = (*sessionListIt)->getFirstChild()->getAttribute("exterAddr");
        IPAddress routerAddr1 = IPAddress(exterAddr);
        exterAddr = (*sessionListIt)->getLastChild()->getAttribute("exterAddr");
        IPAddress routerAddr2 = IPAddress(exterAddr);
        if (isInIPTable(_rt, routerAddr1) == -1 && isInIPTable(_rt, routerAddr2) == -1)
        {
            continue;
        }
        IPAddress peerAddr;
        if (isInIPTable(_rt, routerAddr1) != -1)
        {
            peerAddr = routerAddr2;
            delayTab[3] += atoi((*sessionListIt)->getAttribute("id"));
        }
        else
        {
            peerAddr = routerAddr1;
            delayTab[3] += atoi((*sessionListIt)->getAttribute("id")) + 0.5;
        }
        if (peerAddr.isUnspecified())
        {
            error("In BGPRouting.xml :No valid external address for session ID : %s",(*sessionListIt)->getAttribute("id"));
        }

        BGP::SessionID newSessionID = createSession(BGP::EGP, peerAddr.str().c_str());
        _BGPSessions[newSessionID]->setTimers(delayTab);
        TCPSocket* socketListenEGP = new TCPSocket();
        _BGPSessions[newSessionID]->setSocketListen(socketListenEGP);
    }
}

std::vector<const char *>  BGPRouting::loadASConfig(cXMLElementList& ASConfig)
{
    //create deny Lists
    std::vector<const char *> routerInSameASList;

    for (cXMLElementList::iterator ASConfigIt = ASConfig.begin(); ASConfigIt != ASConfig.end(); ASConfigIt++)
    {
        std::string nodeName = (*ASConfigIt)->getTagName();
        if (nodeName == "Router")
        {
            if (isInIPTable(_rt, IPAddress((*ASConfigIt)->getAttribute("interAddr"))) == -1)
            {
                routerInSameASList.push_back((*ASConfigIt)->getAttribute("interAddr"));
            }
            continue;
        }
        if (nodeName == "DenyRoute" || nodeName == "DenyRouteIN" || nodeName == "DenyRouteOUT")
        {
            BGP::RoutingTableEntry* entry = new BGP::RoutingTableEntry();
            entry->setDestinationID((*ASConfigIt)->getAttribute("Address"));
            entry->setAddressMask((*ASConfigIt)->getAttribute("Netmask"));
            if (nodeName == "DenyRouteIN")
            {
                _prefixListIN.push_back(entry);
            }
            else if (nodeName == "DenyRouteOUT")
            {
                _prefixListOUT.push_back(entry);
            }
            else
            {
                _prefixListIN.push_back(entry);
                _prefixListOUT.push_back(entry);
            }
        }
        else if (nodeName == "DenyAS" || nodeName == "DenyASIN" || nodeName == "DenyASOUT")
        {
            BGP::ASID ASCur = atoi((*ASConfigIt)->getNodeValue());
            if (nodeName == "DenyASIN")
            {
                _ASListIN.push_back(ASCur);
            }
            else if (nodeName == "DenyASOUT")
            {
                _ASListOUT.push_back(ASCur);
            }
            else
            {
                _ASListIN.push_back(ASCur);
                _ASListOUT.push_back(ASCur);
            }
        }
        else
        {
            error("error in BGPConfig.xml : unknown element named '%s' for AS %u", nodeName.c_str(), _myAS);
        }
    }
    return routerInSameASList;
}

bool BGPRouting::loadConfigFromXML (const char * filename)
{
    cXMLElement* bgpConfig = ev.getXMLDocument (filename);
    if (bgpConfig == 0)
    {
        error ("Cannot read BGP configuration from file: %s", filename);
    }

    // load bgp timer parameters informations
    simtime_t delayTab[BGP::NB_TIMERS];
    cXMLElement* paramNode = bgpConfig->getElementByPath("TimerParams");
    if (paramNode == 0)
    {
        error ("In BGPRouting.xml :No configuration for BGP timer parameters");
    }
    cXMLElementList timerConfig = paramNode->getChildren();
    loadTimerConfig(timerConfig, delayTab);

    //find my AS
    cXMLElementList ASList = bgpConfig->getElementsByTagName("AS");
    unsigned char routerPosition;
    _myAS = findMyAS(ASList, _rt, &routerPosition);
    if (_myAS == 0)
    {
        error("In BGPRouting.xml : No AS configuration for Router ID: %s", _rt->getRouterId().str().c_str());
    }

    // load EGP Session informations
    cXMLElementList sessionList = bgpConfig->getElementsByTagName("Session");
    simtime_t saveStartDelay = delayTab[3];
    loadSessionConfig(sessionList, delayTab);
    delayTab[3] = saveStartDelay;

    // load AS information
    char ASXPath[32];
    sprintf(ASXPath, "AS[@id='%d']", _myAS);

    cXMLElement* ASNode = bgpConfig->getElementByPath(ASXPath);
    std::vector<const char *> routerInSameASList;
    if (ASNode != 0)
    {
        cXMLElementList ASConfig = ASNode->getChildren();
        routerInSameASList = loadASConfig(ASConfig);
    }
    else
    {
        error ("In BGPRouting.xml : No configuration for AS ID: %d", _myAS);
    }

    //create IGP Session(s)
    if (routerInSameASList.size())
    {
        unsigned int routerPeerPosition = 1;
        delayTab[3] += sessionList.size()*2;
        for (std::vector<const char *>::iterator it = routerInSameASList.begin(); it != routerInSameASList.end(); it++, routerPeerPosition++)
        {
            BGP::SessionID newSessionID;
            TCPSocket* socketListenIGP = new TCPSocket();
            newSessionID = createSession(BGP::IGP, (*it));
            delayTab[3] += calculateStartDelay(routerInSameASList.size(), routerPosition, routerPeerPosition);
            _BGPSessions[newSessionID]->setTimers(delayTab);
            _BGPSessions[newSessionID]->setSocketListen(socketListenIGP);
        }
    }

    return true;
}

unsigned int calculateStartDelay(int rtListSize, unsigned char rtPosition, unsigned char rtPeerPosition)
{
    unsigned int startDelay = 0;
    if (rtPeerPosition == 1)
    {
        if (rtPosition == 1)
        {
            startDelay = 1;
        }
        else
        {
            startDelay = (rtPosition-1)*2;
        }
        return startDelay;
    }

    if (rtPosition < rtPeerPosition)
    {
        startDelay = 2;
    }
    else if (rtPosition > rtPeerPosition)
    {
        startDelay = (rtListSize-1)*2 - 2*(rtPeerPosition-2);
    }
    else
    {
        startDelay = (rtListSize-1)*2 +1;
    }
    return startDelay;
}

BGP::SessionID BGPRouting::createSession(BGP::type typeSession, const char* peerAddr)
{
    BGPSession*         newSession = new BGPSession(*this);
    BGP::SessionID      newSessionId;
    BGP::SessionInfo    info;

    info.sessionType = typeSession;
    info.ASValue = _myAS;
    info.routerID = _rt->getRouterId();
    info.peerAddr.set(peerAddr);
    if (typeSession == BGP::EGP)
    {
        info.linkIntf = _rt->getInterfaceForDestAddr(peerAddr);
        if (info.linkIntf == 0)
        {
            error("BGPRouting.xml : No configuration interface for peer address: %s", peerAddr);
        }
        info.sessionID = info.peerAddr.getInt() + info.linkIntf->ipv4Data()->getIPAddress().getInt();
    }
    else
    {
        info.sessionID = info.peerAddr.getInt() + info.routerID.getInt();
    }
    newSessionId = info.sessionID;
    newSession->setInfo(info);
    _BGPSessions[newSessionId] = newSession;

    return newSessionId;
}


BGP::SessionID  findIdFromPeerAddr(std::map<BGP::SessionID, BGPSession*> sessions, IPAddress peerAddr)
{
    for (std::map<BGP::SessionID, BGPSession*>::iterator sessionIterator = sessions.begin();
        sessionIterator != sessions.end(); sessionIterator ++)
    {
        if ((*sessionIterator).second->getPeerAddr().equals(peerAddr))
        {
            return (*sessionIterator).first;
        }
    }
    return -1;
}

/*delete BGP Routing entry, if the route deleted correctly return true, false else*/
bool BGPRouting::deleteBGPRoutingEntry(BGP::RoutingTableEntry* entry){
    for (std::vector<BGP::RoutingTableEntry*>::iterator it = _BGPRoutingTable.begin();
         it != _BGPRoutingTable.end (); it++)
    {
        if (((*it)->getHost().getInt () & (*it)->getNetmask().getInt ()) ==
            (entry->getHost().getInt () & entry->getNetmask().getInt ()))
        {
            _BGPRoutingTable.erase(it);
            _rt->deleteRoute(entry);
            return true;
        }
    }
    return false;
}

/*return index of the IP table if the route is found, -1 else*/
int isInIPTable(IRoutingTable* rtTable, IPAddress addr)
{
    for (int i = 0; i < rtTable->getNumRoutes(); i++)
    {
        const IPRoute* entry = rtTable->getRoute(i);
        if (entry->getHost().getInt () == addr.getInt())
        {
            return i;
        }
    }
    return -1;
}

BGP::SessionID  findIdFromSocketConnId(std::map<BGP::SessionID, BGPSession*> sessions, int connId)
{
    for (std::map<BGP::SessionID, BGPSession*>::iterator sessionIterator = sessions.begin();
        sessionIterator != sessions.end(); sessionIterator ++)
    {
        TCPSocket* socket = (*sessionIterator).second->getSocket();
        if (socket->getConnectionId() == connId)
        {
            return (*sessionIterator).first;
        }
    }
    return -1;
}

/*return index of the table if the route is found, -1 else*/
unsigned long BGPRouting::isInTable(std::vector<BGP::RoutingTableEntry*> rtTable, BGP::RoutingTableEntry* entry)
{
    for (unsigned long i = 0; i <rtTable.size(); i++)
    {
        BGP::RoutingTableEntry* entryCur = rtTable[i];
        if ((entry->getHost().getInt () & entry->getNetmask().getInt ()) ==
            (entryCur->getHost().getInt () & entryCur->getNetmask().getInt ()))
        {
            return i;
        }
    }
    return -1;
}

/*return true if the AS is found, false else*/
bool BGPRouting::isInASList(std::vector<BGP::ASID> ASList, BGP::RoutingTableEntry* entry)
{
    for (std::vector<BGP::ASID>::iterator it = ASList.begin(); it != ASList.end(); it++)
    {
        for (unsigned int i = 0; i < entry->getASCount(); i++)
        {
            if ((*it) == entry->getAS(i))
            {
                return true;
            }
        }
    }
    return false;
}

/*return true if OSPF exists, false else*/
bool OSPFExist(IRoutingTable* rtTable)
{
    for (int i=0; i<rtTable->getNumRoutes(); i++)
    {
        if (rtTable->getRoute(i)->getSource() == IPRoute::OSPF)
        {
            return true;
        }
    }
    return false;
}

unsigned char ASLoopDetection(BGP::RoutingTableEntry* entry, BGP::ASID myAS)
{
    for (unsigned int i=1; i< entry->getASCount(); i++)
    {
        if (myAS == entry->getAS(i))
        {
            return BGP::ASLOOP_DETECTED;
        }
    }
    return BGP::ASLOOP_NO_DETECTED;
}

/*return sessionID if the session is found, -1 else*/
BGP::SessionID BGPRouting::findNextSession(BGP::type type, bool startSession)
{
    BGP::SessionID sessionID = -1;
    for (std::map<BGP::SessionID, BGPSession*>::iterator sessionIterator = _BGPSessions.begin();
        sessionIterator != _BGPSessions.end(); sessionIterator ++)
    {
        if ((*sessionIterator).second->getType() == type && !(*sessionIterator).second->isEstablished())
        {
            sessionID = (*sessionIterator).first;
            break;
        }
    }
    if (startSession == true && type == BGP::IGP && sessionID != -1)
    {
        InterfaceEntry* linkIntf = _rt->getInterfaceForDestAddr(_BGPSessions[sessionID]->getPeerAddr());
        if (linkIntf == 0)
        {
            error("No configuration interface for peer address: %s", _BGPSessions[sessionID]->getPeerAddr().str().c_str());
        }
        _BGPSessions[sessionID]->setlinkIntf(linkIntf);
        _BGPSessions[sessionID]->startConnection();
    }
    return sessionID;
}

