var canvas = document.getElementById('stage');
var ctx = canvas.getContext("2d");
ctx.lineWidth= 2;
 var position = 0;
 setInterval(function (){
   ctx.clearRect(0,0,400,400);
   ctx.fillRect(position,0,20,20);
   position++
   if(position> 200){
      position=0;
   }

 }, 20);
