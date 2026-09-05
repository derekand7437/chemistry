import { rnd, ri, shuffle, modeOf, withDistractors } from "./util.js";

/* ---------- chemistry data ---------- */
var AM={H:1.008,He:4.003,Li:6.94,Be:9.012,B:10.81,C:12.011,N:14.007,O:15.999,F:18.998,Ne:20.180,
Na:22.990,Mg:24.305,Al:26.982,Si:28.085,P:30.974,S:32.06,Cl:35.45,Ar:39.948,K:39.098,Ca:40.078,
Ti:47.867,Cr:51.996,Mn:54.938,Fe:55.845,Co:58.933,Ni:58.693,Cu:63.546,Zn:65.38,Br:79.904,
Ag:107.868,Sn:118.710,I:126.904,Ba:137.327,Pb:207.2};

function parseFormula(f){
  var i=0;
  function num(){ var m=/^\d+/.exec(f.slice(i)); if(!m) return 0; i+=m[0].length; return parseInt(m[0],10); }
  function grp(){
    var out={};
    while(i<f.length){
      var ch=f[i];
      if(ch===")") break;
      if(ch==="("){ i++; var inner=grp(); i++; var n=num()||1;
        for(var k in inner) out[k]=(out[k]||0)+inner[k]*n; continue; }
      var el=/^[A-Z][a-z]?/.exec(f.slice(i));
      if(!el){ i++; continue; }
      i+=el[0].length; var c=num()||1;
      out[el[0]]=(out[el[0]]||0)+c;
    }
    return out;
  }
  return grp();
}
function molarMass(f){ var m=parseFormula(f), t=0; for(var k in m) t+=(AM[k]||0)*m[k]; return t; }
function fmt(f){ return f.replace(/(\d+)/g,"<sub>$1</sub>"); }
function sig(x,n){ if(x===0) return "0"; var d=Math.ceil(Math.log10(Math.abs(x))); var p=Math.max(0,(n||3)-d); return (+x.toFixed(Math.min(p,8))).toString(); }

var CATIONS=[
 {s:"Na",c:1,n:"sodium"},{s:"K",c:1,n:"potassium"},{s:"Li",c:1,n:"lithium"},{s:"Ag",c:1,n:"silver"},
 {s:"Mg",c:2,n:"magnesium"},{s:"Ca",c:2,n:"calcium"},{s:"Ba",c:2,n:"barium"},{s:"Zn",c:2,n:"zinc"},
 {s:"Al",c:3,n:"aluminum"},
 {s:"Fe",c:2,n:"iron(II)",rom:1},{s:"Fe",c:3,n:"iron(III)",rom:1},
 {s:"Cu",c:1,n:"copper(I)",rom:1},{s:"Cu",c:2,n:"copper(II)",rom:1},
 {s:"Pb",c:2,n:"lead(II)",rom:1},{s:"Sn",c:4,n:"tin(IV)",rom:1},
 {s:"Cr",c:3,n:"chromium(III)",rom:1},{s:"Mn",c:2,n:"manganese(II)",rom:1},
 {s:"NH4",c:1,n:"ammonium",poly:1}
];
var ANIONS=[
 {s:"Cl",c:1,n:"chloride"},{s:"Br",c:1,n:"bromide"},{s:"I",c:1,n:"iodide"},{s:"F",c:1,n:"fluoride"},
 {s:"O",c:2,n:"oxide"},{s:"S",c:2,n:"sulfide"},{s:"N",c:3,n:"nitride"},
 {s:"NO3",c:1,n:"nitrate",poly:1},{s:"NO2",c:1,n:"nitrite",poly:1},{s:"OH",c:1,n:"hydroxide",poly:1},
 {s:"SO4",c:2,n:"sulfate",poly:1},{s:"SO3",c:2,n:"sulfite",poly:1},{s:"CO3",c:2,n:"carbonate",poly:1},
 {s:"PO4",c:3,n:"phosphate",poly:1},{s:"HCO3",c:1,n:"hydrogen carbonate",poly:1},
 {s:"ClO3",c:1,n:"chlorate",poly:1},{s:"CN",c:1,n:"cyanide",poly:1},{s:"C2H3O2",c:1,n:"acetate",poly:1},
 {s:"MnO4",c:1,n:"permanganate",poly:1},{s:"CrO4",c:2,n:"chromate",poly:1}
];
function gcd(a,b){ return b? gcd(b,a%b):a; }
function ionicFormula(cat,an){
  var g=gcd(cat.c,an.c), a=an.c/g, b=cat.c/g;
  var p1 = cat.poly && a>1 ? "("+cat.s+")" : cat.s;
  var p2 = an.poly && b>1 ? "("+an.s+")" : an.s;
  return p1+(a>1?a:"")+p2+(b>1?b:"");
}
var COVALENT=[
 ["CO","carbon monoxide"],["CO2","carbon dioxide"],["N2O","dinitrogen monoxide"],
 ["N2O4","dinitrogen tetroxide"],["N2O5","dinitrogen pentoxide"],["NO2","nitrogen dioxide"],
 ["PCl3","phosphorus trichloride"],["PCl5","phosphorus pentachloride"],["SF6","sulfur hexafluoride"],
 ["SO2","sulfur dioxide"],["SO3","sulfur trioxide"],["CCl4","carbon tetrachloride"],
 ["P2O5","diphosphorus pentoxide"],["SiO2","silicon dioxide"],["CS2","carbon disulfide"],
 ["ClF3","chlorine trifluoride"],["N2S3","dinitrogen trisulfide"],["P4O10","tetraphosphorus decoxide"]
];
var EQUATIONS=[
 {t:["C3H8","O2","CO2","H2O"],c:[1,5,3,4],r:2},
 {t:["H2","O2","H2O"],c:[2,1,2],r:2},
 {t:["Fe","O2","Fe2O3"],c:[4,3,2],r:2},
 {t:["CH4","O2","CO2","H2O"],c:[1,2,1,2],r:2},
 {t:["Al","HCl","AlCl3","H2"],c:[2,6,2,3],r:2},
 {t:["N2","H2","NH3"],c:[1,3,2],r:2},
 {t:["KClO3","KCl","O2"],c:[2,2,3],r:1},
 {t:["Na","H2O","NaOH","H2"],c:[2,2,2,1],r:2},
 {t:["C2H6","O2","CO2","H2O"],c:[2,7,4,6],r:2},
 {t:["Mg","N2","Mg3N2"],c:[3,1,1],r:2},
 {t:["CaCO3","CaO","CO2"],c:[1,1,1],r:1},
 {t:["Zn","HCl","ZnCl2","H2"],c:[1,2,1,1],r:2},
 {t:["AgNO3","NaCl","AgCl","NaNO3"],c:[1,1,1,1],r:2},
 {t:["C6H12O6","O2","CO2","H2O"],c:[1,6,6,6],r:2},
 {t:["P4","O2","P4O10"],c:[1,5,1],r:2},
 {t:["NaOH","H2SO4","Na2SO4","H2O"],c:[2,1,1,2],r:2}
];
var STOICH=[
 {eq:[["H2",2],["O2",1]],p:[["H2O",2]]},
 {eq:[["N2",1],["H2",3]],p:[["NH3",2]]},
 {eq:[["CH4",1],["O2",2]],p:[["CO2",1],["H2O",2]]},
 {eq:[["Fe",4],["O2",3]],p:[["Fe2O3",2]]},
 {eq:[["Al",2],["HCl",6]],p:[["AlCl3",2],["H2",3]]},
 {eq:[["C3H8",1],["O2",5]],p:[["CO2",3],["H2O",4]]},
 {eq:[["Zn",1],["HCl",2]],p:[["ZnCl2",1],["H2",1]]},
 {eq:[["KClO3",2]],p:[["KCl",2],["O2",3]]}
];
var MM_SET=["H2O","CO2","NaCl","CaCO3","C6H12O6","NH3","H2SO4","MgO","Al2O3","KNO3","C2H6","Fe2O3","CH4","NaOH","Ca(OH)2","N2O5","CuSO4","O2","C3H8","K2CO3"];
var POLY=[
 ["NH4","1+","ammonium"],["NO3","1\u2212","nitrate"],["NO2","1\u2212","nitrite"],["OH","1\u2212","hydroxide"],
 ["CN","1\u2212","cyanide"],["HCO3","1\u2212","bicarbonate"],["ClO3","1\u2212","chlorate"],["ClO4","1\u2212","perchlorate"],
 ["ClO2","1\u2212","chlorite"],["ClO","1\u2212","hypochlorite"],["MnO4","1\u2212","permanganate"],["C2H3O2","1\u2212","acetate"],
 ["HSO4","1\u2212","bisulfate"],["SO4","2\u2212","sulfate"],["SO3","2\u2212","sulfite"],["CO3","2\u2212","carbonate"],
 ["CrO4","2\u2212","chromate"],["Cr2O7","2\u2212","dichromate"],["O2","2\u2212","peroxide"],["C2O4","2\u2212","oxalate"],
 ["PO4","3\u2212","phosphate"],["PO3","3\u2212","phosphite"],["AsO4","3\u2212","arsenate"],["S2O3","2\u2212","thiosulfate"]
];
var TRENDS=[
 {q:"Which atom has the <b>larger atomic radius</b>?",pairs:[["Na","Cl","Na","Both are in period 3. Sodium has fewer protons pulling on the same shell, so its electron cloud is held less tightly."],["K","Na","K","Potassium is one period lower, so it has an entire extra electron shell."],["Mg","Ba","Ba","Barium is far down group 2, and each period adds a shell."],["O","C","C","Across period 2, oxygen's extra protons pull the same shell in tighter."],["Br","Cl","Br","Bromine sits one period below chlorine."],["Al","P","Al","Both are period 3; phosphorus has the stronger nuclear pull."]]},
 {q:"Which atom has the <b>higher ionization energy</b>?",pairs:[["F","Li","F","Ionization energy climbs up and to the right; fluorine holds its electrons hardest in period 2."],["Ne","Na","Ne","Neon has a full shell, so removing an electron is very expensive. Sodium is happy to lose one."],["N","B","N","Moving right across a period, effective nuclear charge grows."],["He","H","He","Helium's two protons hold a full first shell."],["Cl","S","Cl","Chlorine is one step further right in period 3."],["C","Si","C","Carbon's outer electrons are closer to the nucleus and less shielded."]]},
 {q:"Which atom is <b>more electronegative</b>?",pairs:[["F","O","F","Fluorine is the most electronegative element at 4.0."],["O","S","O","Electronegativity increases going up a group."],["N","C","N","Nitrogen is further right in period 2."],["Cl","Br","Cl","Chlorine sits above bromine in group 17."],["O","Na","O","Nonmetals on the right pull shared electrons; metals on the left release them."],["S","P","S","Sulfur is one step further right in period 3."]]}
];

/* ---------- answer checking ---------- */
function norm(s){ return (s||"").toLowerCase().replace(/[\s‐-―_]/g,""); }
function nearly(a,b){ if(!isFinite(a)||!isFinite(b)) return false; if(b===0) return Math.abs(a)<1e-9; return Math.abs(a-b)/Math.abs(b) <= 0.015; }
function parseNum(s){
  s=(s||"").trim().replace(/,/g,"").replace(/\s/g,"").replace(/[xX×]10\^?/,"e").replace(/e\+/,"e");
  var v=parseFloat(s); return isNaN(v)? NaN : v;
}

/* ---------- generators ---------- */
const GENS={};
var ELEMS=[
 ["H","hydrogen",1,0,1],["He","helium",18,0,null],["Li","lithium",1,1,1],["Be","beryllium",2,1,2],
 ["B","boron",13,0,null],["C","carbon",14,0,null],["N","nitrogen",15,0,-3],["O","oxygen",16,0,-2],
 ["F","fluorine",17,0,-1],["Ne","neon",18,0,null],["Na","sodium",1,1,1],["Mg","magnesium",2,1,2],
 ["Al","aluminum",13,1,3],["Si","silicon",14,0,null],["P","phosphorus",15,0,-3],["S","sulfur",16,0,-2],
 ["Cl","chlorine",17,0,-1],["Ar","argon",18,0,null],["K","potassium",1,1,1],["Ca","calcium",2,1,2],
 ["Fe","iron",8,1,null],["Cu","copper",11,1,null],["Zn","zinc",12,1,2],["Ag","silver",11,1,1],
 ["Br","bromine",17,0,-1],["I","iodine",17,0,-1],["Ba","barium",2,1,2],["Pb","lead",14,1,null],
 ["Sn","tin",14,1,null],["Au","gold",11,1,null],["Ni","nickel",10,1,null],["Mn","manganese",7,1,null]
];
var CAT_SIMPLE=CATIONS.filter(function(c){return !c.rom && !c.poly});
var CAT_ROMAN=CATIONS.filter(function(c){return c.rom});
var AN_SIMPLE=ANIONS.filter(function(a){return !a.poly});
var AN_POLY=ANIONS.filter(function(a){return a.poly});
var POLY8=[["NO3","1\u2212","nitrate"],["SO4","2\u2212","sulfate"],["CO3","2\u2212","carbonate"],["PO4","3\u2212","phosphate"],
 ["OH","1\u2212","hydroxide"],["NH4","1+","ammonium"],["C2H3O2","1\u2212","acetate"],["HCO3","1\u2212","bicarbonate"]];
var COUNTABLE=["H2O","CO2","Ca(NO3)2","Al2(SO4)3","NaOH","C6H12O6","Mg(OH)2","(NH4)2SO4","K3PO4","Fe2O3","H2SO4","Cu(NO3)2","C3H8","NH3"];
var GROUPRULE={1:"1+",2:"2+",13:"3+",15:"3\u2212",16:"2\u2212",17:"1\u2212"};

GENS.atoms=function(opt){
  opt=opt||{};
  var mode=modeOf(opt,["sym2name","name2sym","metal","charge","count","type"]);
  if(mode==="sym2name"){
    var e=rnd(ELEMS);
    return {kind:"mc", prompt:"Which element has the symbol <span class=\"big\">"+e[0]+"</span>?",
      choices:withDistractors(e[1],ELEMS.map(function(x){return x[1]}).filter(function(n){return n!==e[1]}),4),
      correct:e[1], answer:e[1],
      steps:["The symbol "+e[0]+" stands for "+e[1]+".",
             "A symbol is one capital letter, sometimes plus one lowercase letter. CO is carbon and oxygen together; Co is the single element cobalt."]};
  }
  if(mode==="name2sym"){
    var e2=rnd(ELEMS);
    if(opt.mc) return {kind:"mc", prompt:"What is the symbol for <b>"+e2[1]+"</b>?",
      choices:withDistractors(e2[0],ELEMS.map(function(x){return x[0]}).filter(function(n){return n!==e2[0]}),4),
      correct:e2[0], answer:e2[0], steps:[e2[1]+" is written as "+e2[0]+".","Some symbols come from Latin names — sodium is Na from <i>natrium</i>, potassium is K from <i>kalium</i>."]};
    return {kind:"text", prompt:"Write the symbol for <b>"+e2[1]+"</b>.", mono:true, placeholder:"e.g. Na",
      check:function(v){ return v.replace(/\s/g,"")===e2[0]; }, near:function(v){ return norm(v)===norm(e2[0]); }, answer:e2[0],
      steps:[e2[1]+" is "+e2[0]+".","Capitalization is part of the answer: the first letter is always capital, the second never is."]};
  }
  if(mode==="metal"){
    var e3=rnd(ELEMS);
    return {kind:"mc", prompt:"Is <b>"+e3[1]+"</b> a metal or a nonmetal?", choices:["metal","nonmetal"],
      correct:e3[3]?"metal":"nonmetal", answer:e3[3]?"metal":"nonmetal",
      steps:[e3[1]+" ("+e3[0]+") sits in group "+e3[2]+", which is on the "+(e3[3]?"left or middle":"right")+" of the table.",
             "Metals fill the left and middle. Nonmetals are the small block in the upper right, past the staircase."]};
  }
  if(mode==="charge"){
    var pool=ELEMS.filter(function(x){return x[4]!==null && GROUPRULE[x[2]]});
    var e4=rnd(pool), ans=GROUPRULE[e4[2]];
    return {kind:"mc", prompt:"What charge does an ion of <b>"+e4[1]+"</b> have?",
      choices:withDistractors(ans,["1+","2+","3+","1\u2212","2\u2212","3\u2212"].filter(function(c){return c!==ans}),4),
      correct:ans, answer:ans,
      steps:[e4[1]+" is in group "+e4[2]+".",
             e4[3]? "Metals lose electrons, so the charge is positive — group 1 gives 1+, group 2 gives 2+, group 13 gives 3+."
                  : "Nonmetals gain electrons, so the charge is negative — group 17 gives 1\u2212, group 16 gives 2\u2212, group 15 gives 3\u2212.",
             "So "+e4[1]+" forms "+e4[0]+"<sup>"+ans+"</sup>."]};
  }
  if(mode==="type"){
    var ionic=Math.random()<0.5, f, why;
    if(ionic){ var c=rnd(CAT_SIMPLE), a=rnd(AN_SIMPLE); f=ionicFormula(c,a); why="It starts with "+c.n+", a metal, joined to a nonmetal. Metal + nonmetal is always ionic — the metal hands its electrons over."; }
    else { f=rnd(COVALENT)[0]; why="Every element in it is a nonmetal, and two nonmetals share electrons instead of transferring them."; }
    return {kind:"mc", prompt:"Is <span class=\"big\">"+fmt(f)+"</span> ionic or covalent?", choices:["ionic","covalent"],
      correct:ionic?"ionic":"covalent", answer:ionic?"ionic":"covalent",
      steps:[why,"This is the first question to ask about any compound, because the two types are named by completely different rules."]};
  }
  var cf=rnd(COUNTABLE), map=parseFormula(cf), keys=Object.keys(map), el=rnd(keys), n=map[el];
  var elName=(ELEMS.filter(function(x){return x[0]===el})[0]||[el,el])[1];
  return {kind:"text", prompt:"How many <b>"+elName+"</b> atoms are in one unit of <span class=\"big\">"+fmt(cf)+"</span>?",
    unit:"atoms", placeholder:"a number",
    check:function(v){ return parseNum(v)===n; }, answer:String(n),
    steps:["A subscript counts only the symbol right in front of it, and no subscript means one.",
           /\(/.test(cf)? "A subscript outside parentheses multiplies everything inside them." : "Read the formula left to right and tally each symbol.",
           "There are "+n+" "+elName+" atom"+(n===1?"":"s")+" in "+cf+"."]};
};

GENS.polyions=function(opt){
  opt=opt||{};
  var ion=rnd(POLY8), toName=modeOf(opt,["f2n","n2f"])==="f2n";
  if(toName) return {kind:"mc", prompt:"What is the name of <span class=\"big\">"+fmt(ion[0])+"<sup>"+ion[1]+"</sup></span>?",
    choices:withDistractors(ion[2],POLY8.map(function(x){return x[2]}).filter(function(n){return n!==ion[2]}),4),
    correct:ion[2], answer:ion[2],
    steps:["That ion is "+ion[2]+", charge "+ion[1]+".","Polyatomic ions travel as one package — the whole group carries a single charge."]};
  return {kind:"mc", prompt:"Which formula is <b>"+ion[2]+"</b>?",
    choices:withDistractors(fmt(ion[0])+"<sup>"+ion[1]+"</sup>",POLY8.filter(function(x){return x[2]!==ion[2]}).map(function(x){return fmt(x[0])+"<sup>"+x[1]+"</sup>"}),4),
    correct:fmt(ion[0])+"<sup>"+ion[1]+"</sup>", answer:ion[0]+" ("+ion[1]+")",
    steps:[ion[2]+" is "+ion[0]+" with a charge of "+ion[1]+".","These have to be memorized — there is no rule that produces them."]};
};

GENS.naming=function(opt){
  opt=opt||{};
  var mode=modeOf(opt,["i2n","n2i","c2n","n2c"]);
  if(mode==="i2n"||mode==="n2i"){
    var cat = opt.roman? rnd(CAT_ROMAN) : (opt.simple||opt.poly)? rnd(CAT_SIMPLE) : rnd(CATIONS);
    var an  = opt.poly? rnd(AN_POLY) : opt.simple? rnd(AN_SIMPLE) : rnd(ANIONS);
    var f=ionicFormula(cat,an), name=cat.n+" "+an.n;
    if(mode==="i2n"){
      var q={kind:"text", prompt:"Name this ionic compound.", big:fmt(f), placeholder:"e.g. sodium chloride",
        check:function(v){ return norm(v)===norm(name); }, answer:name,
        steps:["The metal keeps its name: "+cat.n+(cat.rom?", and the Roman numeral is the charge it must have here.":"."),
               "The nonmetal part is "+an.n+(an.poly?" — a polyatomic ion, which keeps its own name.":" — a single element, so its ending changes to -ide."),
               "Put them in that order: "+name+"."]};
      if(opt.mc){
        var wrong=[];
        for(var i=0;i<12 && wrong.length<3;i++){
          var w=rnd(CATIONS).n+" "+rnd(ANIONS).n;
          if(w!==name && wrong.indexOf(w)<0) wrong.push(w);
        }
        q={kind:"mc", prompt:"What is the name of this compound?", big:fmt(f), choices:shuffle([name].concat(wrong)),
           correct:name, answer:name, steps:q.steps};
      }
      return q;
    }
    return {kind:"text", prompt:"Write the formula for this compound.", big:name, placeholder:"e.g. NaCl", mono:true,
      check:function(v){ return v.replace(/\s/g,"")===f; }, near:function(v){ return norm(v)===norm(f); }, answer:f,
      steps:[cat.n+" is "+cat.s+"<sup>"+cat.c+"+</sup> and "+an.n+" is "+an.s+"<sup>"+an.c+"\u2212</sup>.",
             "Criss-cross: each charge becomes the other ion's subscript, then reduce if they share a factor. The compound has to end up neutral.",
             "Formula: "+f+(/\(/.test(f)?" — parentheses because the polyatomic ion needs a subscript above 1.":".")]};
  }
  var pair=rnd(COVALENT);
  if(mode==="c2n"){
    var qq={kind:"text", prompt:"Name this covalent compound.", big:fmt(pair[0]), placeholder:"e.g. carbon dioxide",
      check:function(v){ return norm(v)===norm(pair[1]); }, answer:pair[1],
      steps:["Two nonmetals, so count the atoms and use prefixes.",
             "Skip mono- on the first element only, and drop a prefix's last vowel before oxide.",
             "Name: "+pair[1]+"."]};
    if(opt.mc) qq={kind:"mc", prompt:"What is the name of this compound?", big:fmt(pair[0]),
      choices:withDistractors(pair[1],COVALENT.map(function(x){return x[1]}).filter(function(n){return n!==pair[1]}),4),
      correct:pair[1], answer:pair[1], steps:qq.steps};
    return qq;
  }
  return {kind:"text", prompt:"Write the formula for this covalent compound.", big:pair[1], placeholder:"e.g. CO2", mono:true,
    check:function(v){ return v.replace(/\s/g,"")===pair[0]; }, near:function(v){ return norm(v)===norm(pair[0]); }, answer:pair[0],
    steps:["Each prefix is a subscript: mono 1, di 2, tri 3, tetra 4, penta 5, hexa 6.",
           "No prefix on the first element means one atom of it.",
           "Formula: "+pair[0]+"."]};
};

GENS.balancing=function(opt){
  opt=opt||{};
  var pool=EQUATIONS;
  if(opt.easy) pool=EQUATIONS.filter(function(e){ return e.t.length<=3 && Math.max.apply(null,e.c)<=4; });
  if(opt.combustion) pool=EQUATIONS.filter(function(e){ return e.t.length===4 && e.t[2]==="CO2"; });
  var e=rnd(pool);
  return {kind:"coefs", prompt:"Balance the equation. Enter every coefficient, including 1.", terms:e.t, react:e.r,
    answer:e.c.join(", "),
    check:function(vals){ for(var i=0;i<e.c.length;i++){ if(parseInt(vals[i],10)!==e.c[i]) return false; } return true; },
    steps:["Count each element on both sides before you change anything.",
           opt.combustion? "Balance carbon first, then hydrogen, and let oxygen fall out last — it appears in both products."
                         : "Balance the element that appears in the fewest places first, and save oxygen for last.",
           "Coefficients: "+e.c.join(", ")+"."]};
};

GENS.moles=function(opt){
  opt=opt||{};
  var set = opt.simple? ["H2O","CO2","NaCl","MgO","CH4","NH3","O2","KNO3","CaCO3"] : MM_SET;
  var f=rnd(set), mm=molarMass(f);
  var mode=modeOf(opt,["mm","g2mol","mol2g","mol2part"]);
  if(mode==="mm") return {kind:"text", prompt:"What is the molar mass of this compound?", big:fmt(f), unit:"g/mol", placeholder:"g/mol",
    check:function(v){ return nearly(parseNum(v),mm); }, answer:sig(mm,5)+" g/mol",
    steps:["Look up each element's atomic mass and multiply it by that element's subscript.",
           "Add every piece — atoms inside parentheses get multiplied through first.",
           "Molar mass = "+sig(mm,5)+" g/mol."]};
  if(mode==="g2mol"){
    var g=ri(5,250)+Math.round(Math.random()*10)/10, mol=g/mm;
    return {kind:"text", prompt:"How many moles are in "+g+" g of "+f+"?", big:g+" g "+fmt(f), unit:"mol", placeholder:"mol",
      check:function(v){ return nearly(parseNum(v),mol); }, answer:sig(mol,4)+" mol",
      steps:["First the molar mass of "+f+": "+sig(mm,5)+" g/mol.",
             "Grams to moles means dividing: "+g+" ÷ "+sig(mm,5)+".",
             "= "+sig(mol,4)+" mol."]};
  }
  if(mode==="mol2g"){
    var mol2=ri(5,60)/10, g2=mol2*mm;
    return {kind:"text", prompt:"What is the mass of "+mol2+" mol of "+f+"?", big:mol2+" mol "+fmt(f), unit:"g", placeholder:"grams",
      check:function(v){ return nearly(parseNum(v),g2); }, answer:sig(g2,4)+" g",
      steps:["Molar mass of "+f+" = "+sig(mm,5)+" g/mol.",
             "Moles to grams means multiplying: "+mol2+" × "+sig(mm,5)+".",
             "= "+sig(g2,4)+" g."]};
  }
  if(mode==="stp"){
    var molS=ri(5,60)/10, L=molS*22.4;
    return {kind:"text", prompt:"What volume does "+molS+" mol of gas fill at STP?", big:"1 mol = 22.4 L at STP", unit:"L", placeholder:"litres",
      check:function(v){ return nearly(parseNum(v),L); }, answer:sig(L,4)+" L",
      steps:["At STP — 0 °C and 1 atm — one mole of any gas fills 22.4 L.",
             molS+" × 22.4 = "+sig(L,4)+" L."]};
  }
  var mol3=ri(5,40)/10, part=mol3*6.022e23;
  return {kind:"text", prompt:"How many molecules are in "+mol3+" mol of "+f+"?", big:mol3+" mol "+fmt(f), unit:"molecules", placeholder:"e.g. 1.2e24",
    check:function(v){ return nearly(parseNum(v),part); }, answer:part.toExponential(3)+" molecules",
    steps:["One mole is 6.022 × 10²³ of anything.",
           mol3+" × 6.022 × 10²³ = "+part.toExponential(3)+".",
           "Type it as "+part.toExponential(3).replace("e+","e")+" — the box reads e-notation."]};
};

GENS.stoich=function(opt){
  opt=opt||{};
  var r=rnd(STOICH), A=rnd(r.eq), B=rnd(r.p);
  var eqHtml=r.eq.map(function(x){return (x[1]>1?x[1]+" ":"")+fmt(x[0])}).join(" + ")+" → "+r.p.map(function(x){return (x[1]>1?x[1]+" ":"")+fmt(x[0])}).join(" + ");
  var mode=modeOf(opt,["mol2mol","g2g"]);
  if(mode==="mol2mol"){
    var molA=ri(10,90)/10, molB=molA*B[1]/A[1];
    return {kind:"text", prompt:"How many moles of "+B[0]+" form from "+molA+" mol of "+A[0]+"?", big:eqHtml, unit:"mol "+B[0], placeholder:"mol",
      check:function(v){ return nearly(parseNum(v),molB); }, answer:sig(molB,4)+" mol",
      steps:["The coefficients give the ratio: "+B[1]+" "+B[0]+" for every "+A[1]+" "+A[0]+".",
             molA+" × ("+B[1]+"/"+A[1]+")",
             "= "+sig(molB,4)+" mol "+B[0]+"."]};
  }
  var mmA=molarMass(A[0]), mmB=molarMass(B[0]), gA=ri(10,200), gB=(gA/mmA)*(B[1]/A[1])*mmB;
  return {kind:"text", prompt:"How many grams of "+B[0]+" form from "+gA+" g of "+A[0]+"?", big:eqHtml, unit:"g "+B[0], placeholder:"grams",
    check:function(v){ return nearly(parseNum(v),gB); }, answer:sig(gB,4)+" g",
    steps:["Grams → moles: "+gA+" ÷ "+sig(mmA,5)+" = "+sig(gA/mmA,4)+" mol "+A[0]+".",
           "Mole ratio "+B[1]+"/"+A[1]+": "+sig((gA/mmA)*(B[1]/A[1]),4)+" mol "+B[0]+".",
           "Moles → grams: × "+sig(mmB,5)+" = "+sig(gB,4)+" g."]};
};

GENS.gases=function(opt){
  opt=opt||{};
  var mode=modeOf(opt,["boyle","charles","gaylussac","ideal","stp"]);
  if(mode==="boyle"){
    var p1=ri(1,4), v1=ri(2,12), p2=ri(1,6); while(p2===p1) p2=ri(1,6);
    var v2=p1*v1/p2;
    return {kind:"text", prompt:"A gas at "+p1+" atm fills "+v1+" L. At the same temperature the pressure changes to "+p2+" atm. What is the new volume?",
      big:"P₁V₁ = P₂V₂", unit:"L", placeholder:"litres",
      check:function(v){ return nearly(parseNum(v),v2); }, answer:sig(v2,4)+" L",
      steps:["Temperature is fixed and only P and V change — Boyle's law.",
             "V₂ = P₁V₁/P₂ = ("+p1+" × "+v1+") / "+p2+".",
             "= "+sig(v2,4)+" L, "+(v2>v1?"bigger, which is right because the pressure dropped.":"smaller, which is right because the pressure rose.")]};
  }
  if(mode==="charles"){
    var vv1=ri(2,10), t1=ri(0,80), t2=ri(90,300), K1=t1+273.15, K2=t2+273.15, vv2=vv1*K2/K1;
    return {kind:"text", prompt:"A balloon holds "+vv1+" L at "+t1+" °C. It is heated to "+t2+" °C at constant pressure. What is the new volume?",
      big:"V₁/T₁ = V₂/T₂", unit:"L", placeholder:"litres",
      check:function(v){ return nearly(parseNum(v),vv2); }, answer:sig(vv2,4)+" L",
      steps:["Kelvin first, always: "+sig(K1,5)+" K and "+sig(K2,5)+" K.",
             "V₂ = V₁ × T₂/T₁ = "+vv1+" × "+sig(K2,5)+"/"+sig(K1,5)+".",
             "= "+sig(vv2,4)+" L — bigger, because heating expands a gas."]};
  }
  if(mode==="gaylussac"){
    var pp1=ri(1,5), tt1=ri(10,60), tt2=ri(120,400), KK1=tt1+273.15, KK2=tt2+273.15, pp2=pp1*KK2/KK1;
    return {kind:"text", prompt:"A sealed steel cylinder reads "+pp1+" atm at "+tt1+" °C. It is heated to "+tt2+" °C. What is the new pressure?",
      big:"P₁/T₁ = P₂/T₂", unit:"atm", placeholder:"atm",
      check:function(v){ return nearly(parseNum(v),pp2); }, answer:sig(pp2,4)+" atm",
      steps:["Sealed and rigid means the volume cannot change — Gay-Lussac's law.",
             "Kelvin: "+sig(KK1,5)+" K and "+sig(KK2,5)+" K, then P₂ = P₁ × T₂/T₁.",
             "= "+sig(pp2,4)+" atm."]};
  }
  if(mode==="ideal"){
    var n=ri(1,6)/2, T=ri(250,500), V=ri(5,40), P=n*0.0821*T/V;
    return {kind:"text", prompt:"What pressure does "+n+" mol of gas exert in a "+V+" L container at "+T+" K?",
      big:"PV = nRT", unit:"atm", placeholder:"atm",
      check:function(v){ return nearly(parseNum(v),P); }, answer:sig(P,4)+" atm",
      steps:["Rearrange to P = nRT/V.",
             "P = ("+n+" × 0.0821 × "+T+") / "+V+".",
             "= "+sig(P,4)+" atm. R is 0.0821 when pressure is in atm and volume in litres."]};
  }
  var molSTP=ri(5,60)/10, L2=molSTP*22.4;
  return {kind:"text", prompt:"What volume does "+molSTP+" mol of gas occupy at STP?", big:"1 mol = 22.4 L at STP", unit:"L", placeholder:"litres",
    check:function(v){ return nearly(parseNum(v),L2); }, answer:sig(L2,4)+" L",
    steps:["STP is 273 K and 1 atm, where a mole of any gas fills 22.4 L.",
           molSTP+" × 22.4 = "+sig(L2,4)+" L."]};
};

GENS.trends=function(opt){
  opt=opt||{};
  var t = (opt.which!=null)? TRENDS[opt.which] : rnd(TRENDS);
  var pr=rnd(t.pairs), opts=[pr[0],pr[1]];
  if(Math.random()<0.5) opts.reverse();
  return {kind:"mc", prompt:t.q, choices:opts, correct:pr[2], answer:pr[2],
    steps:[pr[3],"Radius grows down and to the left. Ionization energy and electronegativity grow up and to the right."]};
};

/* ---------- the 21-day path ---------- */
var PLAN=[
{t:"Atoms and their symbols", min:"8 min", goal:"Read any chemical symbol and say which element it is.",
 teach:[["Everything is made of atoms","Every object around you is built out of about a hundred basic ingredients called <b>elements</b>. An atom is the smallest piece of an element that is still that element."],
 ["The periodic table is the ingredient list","Each box holds one element: its symbol, its name, and its mass. The rows are called periods and the columns are called groups &mdash; elements in the same column behave alike."],
 ["Symbols follow a strict format","One capital letter, sometimes followed by one lowercase letter. That rule carries meaning: <b>CO</b> is carbon joined to oxygen, but <b>Co</b> is the single element cobalt."]],
 tip:"A few symbols come from Latin names instead of English ones: Na sodium, K potassium, Fe iron, Cu copper, Ag silver, Au gold, Pb lead, Sn tin.",
 drill:{topic:"atoms", opt:{modes:["sym2name","name2sym"], mc:true}, target:6}},

{t:"Metals, nonmetals, and ions", min:"9 min", goal:"Predict the charge an atom takes just from where it sits on the table.",
 teach:[["The table has two halves","Metals fill the left and middle &mdash; shiny, bendable, they conduct electricity. Nonmetals are the smaller block in the upper right, past the staircase."],
 ["Atoms trade electrons to get comfortable","An atom is stable when its outer shell is full, so it gains or loses a few electrons to get there. Losing negative electrons leaves a <b>positive</b> ion; gaining them makes a <b>negative</b> ion."],
 ["The column tells you the charge","Group 1 forms 1+, group 2 forms 2+, group 13 forms 3+. Coming from the other side: group 17 forms 1&minus;, group 16 forms 2&minus;, group 15 forms 3&minus;. Group 18 forms nothing &mdash; it is already full."]],
 tip:"Short version: metals lose and turn positive, nonmetals gain and turn negative.",
 drill:{topic:"atoms", opt:{modes:["metal","charge"], mc:true}, target:6}},

{t:"Ionic or covalent?", min:"7 min", goal:"Sort any compound into one of the two families before naming it.",
 teach:[["Metal + nonmetal is ionic","The metal hands its electrons to the nonmetal. Now one is positive and one is negative, and opposite charges hold each other &mdash; that attraction is the bond."],
 ["Nonmetal + nonmetal is covalent","Neither one will give up electrons, so they share a pair instead. Shared pairs make molecules."],
 ["Ask this first, every time","The two families are named by completely different rules. Deciding which family you are in is the real first step of every naming question for the rest of the year."]],
 tip:"Check the first element. If it is a metal, you are in the ionic world.",
 drill:{topic:"atoms", opt:{modes:["type"], mc:true}, target:5}},

{t:"Naming simple ionic compounds", min:"9 min", goal:"Turn a formula like MgCl₂ into a name.",
 teach:[["The metal keeps its name","Nothing changes. NaCl starts with sodium, CaO starts with calcium."],
 ["The nonmetal ends in -ide","Chlorine becomes chlor<b>ide</b>, oxygen becomes ox<b>ide</b>, sulfur becomes sulf<b>ide</b>, nitrogen becomes nitr<b>ide</b>."],
 ["Ignore the subscripts","Ionic names never use prefixes. MgCl₂ is magnesium chloride, not magnesium dichloride &mdash; the subscripts are forced by the charges, so there is no need to say them out loud."]],
 tip:"Two words, always: metal name, then nonmetal name with -ide on the end.",
 drill:{topic:"naming", opt:{modes:["i2n"], simple:true, mc:true}, target:5}},

{t:"Writing ionic formulas", min:"10 min", goal:"Go the other direction: build the formula from the name.",
 teach:[["The compound has to end up neutral","The positives and negatives must cancel exactly. That single requirement decides every subscript."],
 ["Criss-cross the charges","Al³⁺ with O²⁻: the 3 becomes oxygen's subscript and the 2 becomes aluminium's, giving Al₂O₃. Check it: two 3+ is 6+, three 2&minus; is 6&minus;."],
 ["Then reduce","Mg²⁺ with O²⁻ criss-crosses to Mg₂O₂, but both subscripts share a factor of 2, so the real formula is MgO."]],
 tip:"A subscript of 1 is never written. Write NaCl, not Na₁Cl₁.",
 drill:{topic:"naming", opt:{modes:["n2i"], simple:true}, target:5}},

{t:"Polyatomic ions", min:"9 min", goal:"Recognise the eight ions that show up most often.",
 teach:[["Some ions are whole groups of atoms","Nitrate is NO₃⁻ &mdash; one nitrogen and three oxygens carrying a single 1&minus; charge between them. The group travels as one package."],
 ["The eight worth memorising first","nitrate NO₃⁻, sulfate SO₄²⁻, carbonate CO₃²⁻, phosphate PO₄³⁻, hydroxide OH⁻, ammonium NH₄⁺, acetate C₂H₃O₂⁻, bicarbonate HCO₃⁻."],
 ["There is no rule to derive them","This one is pure memory work. Ammonium is the odd one out &mdash; it is the only common positive one."]],
 tip:"Notice the -ate ending on most of them. An -ite ending means the same ion with one less oxygen: sulfate SO₄²⁻, sulfite SO₃²⁻.",
 drill:{topic:"polyions", opt:{}, target:6}},

{t:"Compounds with polyatomic ions", min:"10 min", goal:"Name and build compounds that contain a polyatomic ion.",
 teach:[["The ion keeps its own name","No -ide ending. CaCO₃ is calcium carbonate, NaOH is sodium hydroxide."],
 ["Criss-cross exactly as before","The polyatomic ion's charge behaves like any other charge. Ca²⁺ with NO₃⁻ needs two nitrates."],
 ["Parentheses when you need more than one","Write Ca(NO₃)₂, never CaNO₃₂ &mdash; the parentheses show that the 2 multiplies the whole group."]],
 tip:"Only add parentheses when the subscript on the group is 2 or more. NaOH needs none.",
 drill:{topic:"naming", opt:{modes:["i2n","n2i"], poly:true}, target:5}},

{t:"Transition metals and Roman numerals", min:"9 min", goal:"Work out a metal's charge from the rest of the formula.",
 teach:[["Some metals have more than one charge","Iron can be 2+ or 3+; copper can be 1+ or 2+. The name has to say which one, and it says it with a Roman numeral."],
 ["The numeral is the charge, not the count","Iron(III) chloride means each iron carries a 3+ charge. It does not mean there are three of anything."],
 ["Work backwards from the other ion","In FeCl₃ there are three Cl⁻, so 3&minus; in total. One iron must balance that alone, so iron is 3+ &mdash; iron(III) chloride."]],
 tip:"Numerals I to VII cover everything you will meet: I, II, III, IV, V, VI, VII.",
 drill:{topic:"naming", opt:{modes:["i2n"], roman:true}, target:5}},

{t:"Covalent names use prefixes", min:"9 min", goal:"Name compounds made of two nonmetals.",
 teach:[["Here the numbers do get said","Covalent compounds have no charges to force the subscripts, so the name has to state every count: mono 1, di 2, tri 3, tetra 4, penta 5, hexa 6."],
 ["Skip mono- on the first element only","CO is carbon monoxide, not monocarbon monoxide. But CO₂ is carbon dioxide, and N₂O is dinitrogen monoxide."],
 ["Drop the last vowel before oxide","Monooxide becomes monoxide, pentaoxide becomes pentoxide. It is only for the sound."]],
 tip:"If you see a prefix in a name, you are looking at a covalent compound.",
 drill:{topic:"naming", opt:{modes:["c2n","n2c"]}, target:5}},

{t:"Review: which system?", min:"10 min", goal:"Handle both naming systems mixed together, the way a test does it.",
 teach:[["One question sorts everything","Is the first element a metal? Ionic &mdash; no prefixes, and a Roman numeral if the metal is a transition metal. Not a metal? Covalent &mdash; prefixes on both elements."],
 ["Ionic checklist","Metal name, then the nonmetal with -ide, or the polyatomic ion with its own name. Subscripts come from balancing the charges."],
 ["Covalent checklist","Prefix + first element, prefix + second element with -ide. Every count is spoken aloud."]],
 tip:"When you are stuck, name the two parts separately first, then join them.",
 drill:{topic:"naming", opt:{}, target:6}},

{t:"Counting atoms in a formula", min:"8 min", goal:"Count every atom in a formula, parentheses included.",
 teach:[["A subscript counts one symbol","In H₂O the 2 belongs to hydrogen alone: two hydrogens, one oxygen. No subscript means one."],
 ["Parentheses multiply everything inside","Ca(NO₃)₂ holds one calcium, two nitrogens, and six oxygens &mdash; three oxygens in the group, doubled."],
 ["Coefficients multiply the whole formula","3 H₂O is six hydrogens and three oxygens. You will need this the moment you start balancing tomorrow."]],
 tip:"Write a tally down the side of the page. Every balancing mistake starts as a counting mistake.",
 drill:{topic:"atoms", opt:{modes:["count"]}, target:5}},

{t:"Balancing: the easy ones", min:"10 min", goal:"Balance a short equation by adjusting coefficients.",
 teach:[["Atoms are never created or destroyed","Whatever goes into a reaction comes out of it. So each element must appear the same number of times on both sides."],
 ["You may only change the big numbers in front","Changing H₂O into H₂O₂ makes a different substance &mdash; hydrogen peroxide. Put a coefficient in front instead."],
 ["Count, adjust, recount","List each element with its count on both sides. Fix the worst mismatch, then recount everything, because a fix changes other counts too."]],
 tip:"Balance the element that appears in the fewest places first. Save oxygen for last.",
 drill:{topic:"balancing", opt:{easy:true}, target:4}},

{t:"Balancing combustion reactions", min:"10 min", goal:"Balance a fuel burning in oxygen.",
 teach:[["Combustion has a fixed shape","Something with carbon and hydrogen plus O₂, producing CO₂ and H₂O every single time."],
 ["Carbon, then hydrogen, then oxygen","Carbon only appears in CO₂ and hydrogen only in H₂O, so both are easy. Once they are set, count the oxygens on the right and fill in O₂ last."],
 ["Odd oxygen? Double everything","If the right side needs an odd number of oxygens, put a 2 in front of the whole fuel and start again &mdash; that turns the odd count even."]],
 tip:"C₃H₈ + 5 O₂ → 3 CO₂ + 4 H₂O. Three carbons, eight hydrogens, then ten oxygens on the right means five O₂.",
 drill:{topic:"balancing", opt:{combustion:true}, target:4}},

{t:"Molar mass", min:"9 min", goal:"Find the mass of one mole of any compound.",
 teach:[["A mole is a counting word","Like a dozen, only much bigger: 6.022 × 10²³ of something. Chemists count in moles because atoms are far too small to count one at a time."],
 ["Molar mass is on the periodic table","The decimal number in each box is that element's mass in grams per mole. Carbon is 12.01 g/mol."],
 ["Add up the whole formula","Multiply each element's mass by its subscript, then add. For H₂O: 2(1.008) + 16.00 = 18.02 g/mol."]],
 tip:"Round to two decimals as you go. The answers here allow about 1% of slack.",
 drill:{topic:"moles", opt:{modes:["mm"], simple:true}, target:5}},

{t:"Grams and moles", min:"10 min", goal:"Convert between a mass on a balance and a number of moles.",
 teach:[["Molar mass is the exchange rate","It says how many grams sit in one mole, so it converts between the two directions."],
 ["Grams to moles: divide","36.0 g of water ÷ 18.02 g/mol = 2.00 mol."],
 ["Moles to grams: multiply","2.00 mol of water × 18.02 g/mol = 36.0 g. If you cannot remember which way, watch the units cancel."]],
 tip:"Write the conversion as a fraction with the unwanted unit on the bottom. If the units cancel, the setup is right.",
 drill:{topic:"moles", opt:{modes:["g2mol","mol2g"], simple:true}, target:5}},

{t:"Particles and gas volume", min:"9 min", goal:"Convert moles into molecules, and into litres of gas.",
 teach:[["Avogadro's number","One mole is 6.022 × 10²³ particles &mdash; atoms, molecules, or formula units. Multiply moles by it to get particles; divide to come back."],
 ["Gases at STP","At 0 °C and 1 atm, one mole of <i>any</i> gas fills 22.4 litres. The identity of the gas does not matter."],
 ["The mole sits in the middle","Grams, particles, and litres never convert directly into each other. Everything passes through moles."]],
 tip:"Type big numbers as 1.2e24 &mdash; the answer box reads that as 1.2 × 10²⁴.",
 drill:{topic:"moles", opt:{modes:["mol2part","stp"], simple:true}, target:5}},

{t:"Mole ratios", min:"9 min", goal:"Use a balanced equation to get from one substance to another.",
 teach:[["A balanced equation is a recipe","2 H₂ + O₂ → 2 H₂O reads: two moles of hydrogen and one of oxygen make two of water."],
 ["The coefficients are the ratio","They are the only link between one substance and another. Nothing else in the problem connects them."],
 ["Write the ratio so units cancel","Going from moles of H₂ to moles of H₂O, put H₂O on top: 5.0 mol H₂ × (2 H₂O / 2 H₂) = 5.0 mol H₂O."]],
 tip:"If the equation is not balanced, balance it before anything else. Every number afterwards depends on it.",
 drill:{topic:"stoich", opt:{modes:["mol2mol"]}, target:5}},

{t:"Grams to grams", min:"11 min", goal:"Run the full three-step conversion, the centrepiece of the course.",
 teach:[["The road never changes","grams A → moles A → moles B → grams B. Three steps, in that order, every time."],
 ["Step one and three use molar mass","Divide by A's molar mass on the way in; multiply by B's molar mass on the way out."],
 ["Step two uses the coefficients","That middle step is the only place the balanced equation appears, and it is the step people forget."]],
 tip:"Grams never convert straight into grams. If your work has no mole in the middle, it is wrong.",
 drill:{topic:"stoich", opt:{modes:["g2g"]}, target:4}},

{t:"Boyle's law and Charles's law", min:"9 min", goal:"Predict what happens to a gas when you squeeze it or heat it.",
 teach:[["Kelvin, always","K = °C + 273.15. Celsius in a gas law will hand you a negative volume, which is impossible."],
 ["Boyle: pressure and volume fight","Hold the temperature still and squeezing a gas into half the space doubles its pressure. P₁V₁ = P₂V₂."],
 ["Charles: temperature and volume agree","Hold the pressure still and heating a gas expands it. V₁/T₁ = V₂/T₂."]],
 tip:"Before trusting a number, check the direction. Squeezed harder but your pressure went down? A ratio is upside down.",
 drill:{topic:"gases", opt:{modes:["boyle","charles"]}, target:5}},

{t:"Gay-Lussac and the ideal gas law", min:"10 min", goal:"Handle sealed containers and problems that involve an amount of gas.",
 teach:[["Gay-Lussac: sealed and rigid","If the container cannot change size, heating it raises the pressure instead. P₁/T₁ = P₂/T₂ &mdash; this is why sealed cans burst in a fire."],
 ["The ideal gas law brings in moles","PV = nRT, with R = 0.0821 L·atm/mol·K. Use it when a problem mentions an amount of gas rather than a before and after."],
 ["Match your units to R","Pressure in atmospheres, volume in litres, temperature in kelvin. Any other unit has to be converted first."]],
 tip:"Two states given? Use a before-and-after law. One state and an amount? Use PV = nRT.",
 drill:{topic:"gases", opt:{modes:["gaylussac","ideal","stp"]}, target:5}},

{t:"Periodic trends", min:"9 min", goal:"Explain the table's patterns instead of memorising them.",
 teach:[["Two forces explain everything","More protons pull the electrons in tighter. More inner shells push the outer ones away and shield them from that pull."],
 ["Atomic radius grows down and left","Going down adds whole shells. Going left removes protons, weakening the pull."],
 ["Ionization energy and electronegativity grow up and right","Both measure how strongly an atom holds electrons, so they run opposite to radius. Fluorine, in the top right corner, is the greediest element."]],
 tip:"Cations are always smaller than their atom, anions always bigger. Losing an electron often removes a whole shell.",
 drill:{topic:"trends", opt:{}, target:6}}
];

export const CHEM = {
  id: "chemistry",
  name: "The Chem Bench",
  storeKey: "chembench",
  generators: GENS,
  plan: PLAN,
  formatTerm: fmt,
  bodyClass: "card-body",
  labelClass: "eyebrow",
  finale: "Twenty-one days, first semester covered. The library below never runs out of problems — and the flashcards are still the fastest review before a test.",
  firstVisit: "Start at Day 1 — it assumes you know nothing at all."
};
export { GENS, PLAN, fmt, POLY };
