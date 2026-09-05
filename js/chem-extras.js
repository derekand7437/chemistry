import { $, $$ } from "./util.js";
import { fmt, POLY } from "./chem-content.js";

export function mountExtras(){
  /* ---------- flashcards ---------- */
  $("#cards").innerHTML = POLY.map(function(p){
    var sym=fmt(p[0])+"<sup>"+p[1]+"</sup>";
    return '<button class="fcard" aria-label="Flashcard: '+p[2]+'"><span class="inner">'+
      '<span class="face front"><span class="sym">'+sym+'</span><span class="tag">tap to flip</span></span>'+
      '<span class="face back"><span class="nm">'+p[2]+'</span><span class="tag">charge '+p[1]+'</span></span>'+
      '</span></button>';
  }).join("");
  $$(".fcard").forEach(function(c){ c.addEventListener("click",function(){ c.classList.toggle("flip"); }); });
}
