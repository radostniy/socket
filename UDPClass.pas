unit UDPClass;

interface

uses
  Classes, Windows, Messages, WinSock, SysUtils;


type
  TNotifyEventRecive = procedure(Sender: TObject; Packet: array of byte) of object;
  TNotifyEventError  = procedure(Sender: TObject; ErrTxt: string) of object;

  TUDP = class({TObject}TThread)

  private
    FPort: Integer;
    FReciveAdr: string;
    FLenP: Integer;
    FOnReciveData : TNotifyEventRecive;
    FSocketError : TNotifyEventError;

  protected
    FSock: TSocket;
    FServerAddr: TSockAddrIn;
    FwData: TWSAData;
    FlUDP: Boolean;
    FPacket: array of Byte;
    th1: cardinal;
    h1: integer;
    FTxt: string;

    procedure Execute; override;

  public
    property Port:Integer write FPort;
    property ReciveAdr:string read FReciveAdr;
    property ReciveLenP:Integer read FLenP write FLenP;

    procedure Terminate;

    procedure SendToIP(fip: AnsiString; port: Integer; const Packet: array of byte);

  published
    property OnReciveData  : TNotifyEventRecive read FOnReciveData write FOnReciveData;
    property OnSocketError : TNotifyEventError  read FSocketError  write FSocketError;
  end;


implementation

procedure TUDP.Terminate;
begin
  FlUDP := False;
  closesocket(FSock);
  WSACleanup;

  inherited Terminate;
end;

procedure TUDP.Execute;
var
  Count: LongInt;
  from: TSockAddrIn;
  FromLen: LongInt;
  Packet: array of Byte;

begin
  FreeOnTerminate := True;

  WSAStartup($0002, FwData);
  FSock := socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP);
  FServerAddr.sin_family := AF_INET;
  FServerAddr.sin_port := htons(FPort);
  bind(FSock, FServerAddr, SizeOf(FServerAddr));

  FlUDP := True;

  while {FlUDP} not Terminated do
  begin
    Sleep(1);

    if Length(Packet)<>FLenP then
      SetLength(Packet, FLenP);

    FromLen := SizeOf(from);
    Count := recvfrom(FSock, Packet[0], FLenP, 0, from, FromLen);
    //Count := recv(FSock, Packet[0], FLenP, 0);

//    FTxt := IntToStr(Count);
//    if Assigned(OnSocketError) then
//        OnSocketError(Self, FTxt);

    if (Terminated) or (not FlUDP) then
      Exit;


    if Count = INVALID_SOCKET then
    begin
      FTxt := 'Произошла ошибка сокета: ' + SysErrorMessage(WSAGetLastError);

      if Assigned(OnSocketError) then
        OnSocketError(Self, FTxt);
    end;

    if Count > 0 then
    begin
      FReciveAdr :=  string(inet_ntoa(from.sin_addr));

      if Assigned(FOnReciveData) then
      begin
        //if Length(FPacket)<>Count then
        //  SetLength(FPacket, Count);

        //CopyMemory(@FPacket[0], @Packet[0], Count);
        FOnReciveData(Self, {FPacket}Packet);
      end;
    end;
  end;
end;

procedure TUDP.SendToIP(fip: AnsiString; port: Integer; const Packet: array of byte);
begin
  FServerAddr.sin_port := htons(port);
  FServerAddr.sin_addr.S_addr := inet_addr(PAnsiChar(fip));
  sendto(FSock, Packet[0], SizeOf(Packet), 0, FServerAddr, SizeOf(FServerAddr));
end;


end.
