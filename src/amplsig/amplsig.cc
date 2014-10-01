/*
 AMPL extension DLL that simulates sending SIGINT (Ctrl-C) from another
 process on Windows. It is required since GenerateConsoleCtrlEvent doesn't
 work reliably with CTRL_C_EVENT across processes.
 The DLL can be loaded from AMPL as follows:
   load signal.dll;
 It connects to a socket at 127.0.0.1:$AMPLSIG_PORT and generates
 CTRL_C_EVENT on every SIGINT command read from the socket input stream.

 Copyright (C) 2014 AMPL Optimization Inc

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization Inc disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Victor Zverovich
 */

#include <cstdio>
#include <cstdlib>

#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <winsock2.h>
#include <ws2tcpip.h>

#pragma comment(lib, "Ws2_32.lib")

using namespace std;

BOOL APIENTRY DllMain(HMODULE, DWORD, LPVOID) {
  return TRUE;
}

class Thread {
 private:
  WSADATA wsaData;
  addrinfo* addr;
  SOCKET socket;
  int debugLevel;

  // Do not implement!
  Thread(const Thread&);
  void operator=(const Thread&);

  void fail(const char* funcName, int code);
  void log(const char* message) {
    if (debugLevel > 1)
      fprintf(stderr, "amplsig: %s\n", message);
  }

  Thread();
  ~Thread();

  void doRun();

 public:
  static DWORD WINAPI run(LPVOID);
};

void Thread::fail(const char* funcName, int code) {
  if (debugLevel > 0)
    fprintf(stderr, "amplsig: %s failed: %d\n", funcName, code);
  throw code;
}

Thread::Thread() : addr(0), socket(INVALID_SOCKET), debugLevel(0) {
  const char* debug = getenv("AMPLSIG_DEBUG");
  if (debug) debugLevel = atoi(debug);
  int result = WSAStartup(MAKEWORD(2, 2), &wsaData);
  if (result != 0)
    fail("WSAStartup", result);
}

Thread::~Thread() {
  if (socket != INVALID_SOCKET)
    closesocket(socket);
  if (addr)
    freeaddrinfo(addr);
  WSACleanup();
}

void Thread::doRun() {
  const char* port = getenv("AMPLSIG_PORT");
  if (!port) {
    log("Port is not specified");
    return;
  }

  addrinfo hints = {0};
  hints.ai_family = AF_UNSPEC;
  hints.ai_socktype = SOCK_STREAM;
  hints.ai_protocol = IPPROTO_TCP;
  int result = getaddrinfo("127.0.0.1", port, &hints, &addr);
  if (result != 0)
    fail("getaddrinfo", result);
  socket = ::socket(addr->ai_family,
    addr->ai_socktype, addr->ai_protocol);
  if (socket == INVALID_SOCKET)
    fail("socket", WSAGetLastError());
  log("Connecting");
  result = connect(socket, addr->ai_addr, (int)addr->ai_addrlen);
  if (result == SOCKET_ERROR)
    fail("connect", WSAGetLastError());
  log("Connected");

  enum { BUFFER_LEN = 512 };
  char buffer[BUFFER_LEN];
  do {
    result = recv(socket, buffer, BUFFER_LEN, 0);
    if (result > 0) {
      if (strncmp(buffer, "SIGINT", 6) == 0 &&
          !GenerateConsoleCtrlEvent(CTRL_C_EVENT, 0)) {
        fail("GenerateConsoleCtrlEvent", GetLastError());
      }
    } else if (result != 0)
      fail("recv", WSAGetLastError());
  } while (result > 0);
}

DWORD WINAPI Thread::run(LPVOID) {
  try {
    Thread().doRun();
  } catch (...) {}
  return 0;
}

extern "C" __declspec(dllexport) void funcadd_ASL(void*) {
  CreateThread(NULL, 0, Thread::run, NULL, 0, NULL);
}
