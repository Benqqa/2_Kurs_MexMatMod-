var canvas;
var ctx = canvas.getContext("2d");
ctx.lineWidth= 2;
 var position = 0;
 setInterval(Draw(), 20);
 function Draw(){
   canvas = document.getElementById('stage');
   ctx = canvas.getContext("2d");
   ctx.clearRect(0,0,400,400);
   ctx.fillRect(position,0,20,20);
   position++
   if(position> 200){
      position=0;
   }
 }
