// Setup server

const http = require('http')
const fs = require('fs')

// Loading the index file . html displayed to the client

const server = http.createServer(function(req, res) {
  fs.readFile('index.html', 'utf-8', function(error, content) {
    res.writeHead(200, {"Content-Type": "text/html; charset=utf-8"})
    res.end(content)
  })
})

// Loading socket.io

const io = require('socket.io').listen(server)
// When a client connects, we note it in the console
io.sockets.on('connection', function (socket) {
  console.log('A client is connected!')
});
server.listen(8084)

// End of setuping server

module.exports = server