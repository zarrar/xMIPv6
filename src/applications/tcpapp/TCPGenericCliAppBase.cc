//
// Copyright 2004 Andras Varga
//
// This library is free software, you can redistribute it and/or modify
// it under  the terms of the GNU Lesser General Public License
// as published by the Free Software Foundation;
// either version 2 of the License, or any later version.
// The library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU Lesser General Public License for more details.
//


#include "TCPGenericCliAppBase.h"

#include "GenericAppMsg_m.h"
#include "IPAddressResolver.h"


simsignal_t TCPGenericCliAppBase::connectSignal = SIMSIGNAL_NULL;
simsignal_t TCPGenericCliAppBase::rcvdPkBytesSignal = SIMSIGNAL_NULL;
simsignal_t TCPGenericCliAppBase::sentPkBytesSignal = SIMSIGNAL_NULL;

void TCPGenericCliAppBase::initialize()
{
    cSimpleModule::initialize();
    numSessions = numBroken = packetsSent = packetsRcvd = bytesSent = bytesRcvd = 0;

    //statistics
    connectSignal = registerSignal("connect");
    rcvdPkBytesSignal = registerSignal("rcvdPkBytes");
    sentPkBytesSignal = registerSignal("sentPkBytes");

    emit(connectSignal, 0L);

    WATCH(numSessions);
    WATCH(numBroken);
    WATCH(packetsSent);
    WATCH(packetsRcvd);
    WATCH(bytesSent);
    WATCH(bytesRcvd);

    // parameters
    const char *address = par("address");
    int port = par("port");
    socket.readDataTransferModePar(*this);
    socket.bind(*address ? IPvXAddress(address) : IPvXAddress(), port);

    socket.setCallbackObject(this);
    socket.setOutputGate(gate("tcpOut"));

    setStatusString("waiting");
}

void TCPGenericCliAppBase::handleMessage(cMessage *msg)
{
    if (msg->isSelfMessage())
        handleTimer(msg);
    else
        socket.processMessage(msg);
}

void TCPGenericCliAppBase::connect()
{
    // we need a new connId if this is not the first connection
    socket.renewSocket();

    // connect
    const char *connectAddress = par("connectAddress");
    int connectPort = par("connectPort");

    EV << "issuing OPEN command\n";
    setStatusString("connecting");

    socket.connect(IPAddressResolver().resolve(connectAddress), connectPort);

    numSessions++;
    emit(connectSignal, 1L);
}

void TCPGenericCliAppBase::close()
{
    setStatusString("closing");
    EV << "issuing CLOSE command\n";
    socket.close();
    emit(connectSignal, 0L);
}

void TCPGenericCliAppBase::sendPacket(int numBytes, int expectedReplyBytes, bool serverClose)
{
    EV << "sending " << numBytes << " bytes, expecting " << expectedReplyBytes << (serverClose ? ", and server should close afterwards\n" : "\n");

    GenericAppMsg *msg = new GenericAppMsg("data");
    msg->setByteLength(numBytes);
    msg->setExpectedReplyLength(expectedReplyBytes);
    msg->setServerClose(serverClose);

    socket.send(msg);

    packetsSent++;
    bytesSent += numBytes;
    emit(sentPkBytesSignal, numBytes);
}

void TCPGenericCliAppBase::setStatusString(const char *s)
{
    if (ev.isGUI()) getDisplayString().setTagArg("t", 0, s);
}

void TCPGenericCliAppBase::socketEstablished(int, void *)
{
    // *redefine* to perform or schedule first sending
    EV << "connected\n";
    setStatusString("connected");
}

void TCPGenericCliAppBase::socketDataArrived(int, void *, cPacket *msg, bool)
{
    // *redefine* to perform or schedule next sending
    packetsRcvd++;
    bytesRcvd += msg->getByteLength();
    emit(rcvdPkBytesSignal, (long)(msg->getByteLength()));

    delete msg;
}

void TCPGenericCliAppBase::socketPeerClosed(int, void *)
{
    // close the connection (if not already closed)
    if (socket.getState()==TCPSocket::PEER_CLOSED)
    {
        EV << "remote TCP closed, closing here as well\n";
        close();
    }
}

void TCPGenericCliAppBase::socketClosed(int, void *)
{
    // *redefine* to start another session etc.
    EV << "connection closed\n";
    setStatusString("closed");
}

void TCPGenericCliAppBase::socketFailure(int, void *, int code)
{
    // subclasses may override this function, and add code try to reconnect after a delay.
    EV << "connection broken\n";
    setStatusString("broken");

    numBroken++;
}

void TCPGenericCliAppBase::finish()
{
    EV << getFullPath() << ": opened " << numSessions << " sessions\n";
    EV << getFullPath() << ": sent " << bytesSent << " bytes in " << packetsSent << " packets\n";
    EV << getFullPath() << ": received " << bytesRcvd << " bytes in " << packetsRcvd << " packets\n";

    recordScalar("number of sessions", numSessions);
}
