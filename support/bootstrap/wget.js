// Poor man's wget for Windows.

var fso = new ActiveXObject("Scripting.FileSystemObject");
var filename = WScript.Arguments(1);
if (fso.FileExists(filename))
  fso.DeleteFile(filename);
var http_req = new ActiveXObject("WinHttp.WinHttpRequest.5.1");
http_req.Open("GET", WScript.Arguments(0), /*async=*/false);
http_req.Send();
stream = new ActiveXObject("ADODB.Stream");
stream.Type = 1;
stream.Open();
stream.Write(http_req.ResponseBody);
stream.SaveToFile(filename);
