var express = require("express");
var exec = require("child_process").exec;
var app = express();

// exec("cmd /k hexo clean");
// exec("cmd /k hexo g");

app.use(express.static('public'));

var server = app.listen(3000, '0.0.0.0', function () {
    var host = server.address().address;
    var port = server.address().port;
    console.log("Server listen at http://localhost:%d", port);
}); 
