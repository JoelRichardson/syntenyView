var mmServer = "www.mousemine.org";

//let mgiUrlBase = "http://www.informatics.jax.org/batch/summary";
//let mmUrlBase  = "http://www.mousemine.org/mousemine/portal.do?class=Gene&externalids=";
var mousemine   = new intermine.Service({root: mmServer+'/mousemine/service'});
var margin = {top: 35, right: 10, bottom: 20, left: 30};
var outerWidth = 800;
var outerHeight = 400;
var width = outerWidth - margin.left - margin.right;
var height = outerHeight - margin.top - margin.bottom;
var zoomer = d3.behavior.zoom(); // .on("zoom", zoom); DISABLE
var yscale = null;

var species = []; // list of species data
var spIndex = 0; // index of species for backbone
var cIndex = 0; // index of species for coloring
var colors=null; // color map

var labels = null;

var aName = null; // A genome name
var bName = null; // B genome name

var aSelected = []; // currently selected A chrs
var bSelected = []; // currently selected B chrs

var allBlocks = []; // the synteny blocks

var minRectHeight = 2; // Min size in px for drawing a block. 0 for no minimum
var magnification = 1.0;
var maxBlockSize = 1000000; // Max size in bp of synteny blocks to draw.
var inflationThreshold = 1.1; // only show blocks whose inflation is >= this value (1 == show everything)

var facets = []; // empty to show all
var facetFuncs = {
    "inversions" : function (blk) { return blk.ori === "-" },
    "translocations" : function (blk) { return blk.aChr !== blk.bChr },
    "deletions" : function (blk) { return blk.aChr === null || blk.bChr === null; },
    "inflation" : function (blk) { return blk.inflation >= inflationThreshold },
    "maxBlockSize" : function (blk) { return blk.bLength <= maxBlockSize }
};

var cwidth = 20;
var dur = 1500;
var bwidth = cwidth/2;

var yViewSize = null;
var svg = d3.select("#display")
	.attr("width", outerWidth)
	.attr("height", outerHeight)
      .append("g")
	.attr("transform", "translate(" + margin.left + "," + margin.top + ")")
      .append("g")
	.call(zoomer)
     .append("g")
	;
var overlay = svg.append('rect').attr('width',width).attr('height',height).attr('class','overlay');
var text =  d3.select("#display")
	.append("g")
	.append("text")
	    .text("")
	    .attr("x", 100)
	    .attr("y", outerHeight-10)
	    .attr("fill", "black")
	;

// ------------------------------
d3.select("#showFacets")
    .selectAll('input[type="checkbox"')
    .on("change", setShowFacets);

function setShowFacets () {
    facets = [];
    var cbs = d3.select("#showFacets")
        .selectAll('input[type="checkbox"')
        .each( function () {
            this.checked && facets.push(this.value) 
        });
    redraw();
}

// ------------------------------
d3.select('[name="maxBlockSize"].labelledInput select')
    .on("change", function () { setMaxBlockSize(this.value) });

function setMaxBlockSize(n) {
    maxBlockSize = +n;
    redraw();
}

// ------------------------------
d3.select('[name="minRectHeight"].labelledInput select')
    .on("change", function () { setMinRectHeight(this.value) });

function setMinRectHeight(n) {
    minRectHeight = +n;
    redraw();
}

// ------------------------------
d3.select('[name="inflation"].labelledInput input[type="text"]')
    .on("change", function () { setInflationThreshold(this.value) });

function setInflationThreshold(t) {
    inflationThreshold = +t;
    redraw();
}

// ---------------------------------------------
// ---------------------------------------------
//
function name2text(n){
    if (n.startsWith('mus_musculus_')) {
        return n.replace('mus_musculus_', '');
    }
    else if (n.startsWith('mus_spretus_')) {
        return n.replace('mus_spretus_', '');
    }
    else
        return n.replace('mus_', 'mus ')
}
d3tsv("./output/strainList.tsv").then(function(strains) {
    let s1 = initOptList("#aGenome", strains, s=>s.strain, s=>name2text(s.strain))
    let s2 = initOptList("#bGenome", strains, s=>s.strain, s=>name2text(s.strain))
    s1.on("change", function () {
        let s2strs = strains.slice(s1.node().selectedIndex);
        initOptList("#bGenome", s2strs, s=>s.strain, s=>name2text(s.strain));
        s2.node().selectedIndex = s2strs.length > 1 ? 1 : 0;
        go();
    });
    s2.on("change", go);
    s1[0][0].selectedIndex = 0;
    s2[0][0].selectedIndex = 3;
    go();
});

// ---------------------------------------------
// Initialize the strain selection lists. 
function initOptList( selector, opts, value, text ) {
    let s = d3.select(selector);
    let os = s.selectAll("option").data(opts);
    let ident = d => d;
    value = value || ident;
    text = text || value;
    os.enter().append("option") ;
    os.exit().remove() ;
    os.attr("value", value)
      .text( text ) ;
    return s;
}

function formatLength(bp) {
    let abp = Math.abs(bp);
    if (abp < 1000)
        return `${bp} bp`
    else if (abp < 1000000) 
        return `${(bp/1000).toFixed(1)} kb`
    else
        return `${(bp/1000000).toFixed(1)} Mb`
}

//
// Promisifies a call to d3.tsv.
// Args:
//   url (string) The url of the json resource
// Returns:
//   a promise that resolves to the json object value, or rejects with an error
function d3tsv(url) {
    return new Promise(function(resolve, reject) {
        d3.tsv(url, function(error, tsv){
            error ? reject({ status: error.status, statusText: error.statusText}) : resolve(tsv);
        })  
    }); 
}

function go () {
    aName = d3.select("#aGenome")[0][0].value;
    bName = d3.select("#bGenome")[0][0].value;
    console.log("GO!", aName, bName);
    let makefn = (a,b) => `./output/${a}-${b}.tsv`;
    Promise.all([
        d3tsv(makefn(aName,aName)),
        d3tsv(makefn(bName,bName)),
        d3tsv(makefn(aName,bName))
    ]).then(function(data) {
        let aChrs = data[0].map(c => [ c.aChr, c.aLength, aName ]);
        let bChrs = data[1].map(c => [ c.bChr, c.bLength, bName ]);
        let abBlks = data[2];
        let allChrs = aChrs.concat(bChrs).sort((a,b) => {
            let ca = a[0];
            let cb = b[0];
            let can = parseInt(ca);
            let cbn = parseInt(cb);
            if (isNaN(can) || isNaN(cbn))
                return ca < cb ? -1 : ca > cb ? 1 : 0;
            else
                return can - cbn;
        });

        let bks = [];
        abBlks.forEach(function(k){
          bks.push({
            name   : ""+k.blockId,
            ori    : k.blockOri,
            aChr   : k.aChr,
            aStart : k.aStart,
            aEnd   : k.aEnd,
            aLength: k.aLength,
            bChr   : k.bChr,
            bStart : k.bStart,
            bEnd   : k.bEnd,
            bLength: k.bLength,
            ids    : k.ids.split(','),
            inflation : 1 / k.blockRatio,
            map    : d3.scale.linear().clamp(true)
          });
        });
        allBlocks = bks;
        aSelected = [];
        bSelected = [];
        species = [];
        sp2chrs = {}
        processChrs( allChrs );
    });
}

function processChrs(data) { 
    // data = Data for each chromosome in order by id
    // (chrs from different species are interleaved)
    //
    // Data = list of: primaryIdentifier, length, organism.shortName
    //
    var sp, chr, len;
    var s2chrs = {};
    var spd, dd; // spd=species data, dd=data for 1 chr
    var sx,sy;
    var data2 = [];
    var maxLen = 0; // for any chromosome of any species

    species = [];
    data.forEach(function(d){
	sp = d[2];      // "species" name
	chr = d[0];     // chromosome name
	len = +d[1];    // chromosome length
        // skip contigs
	if(chr.length > 2) return;
        //
	spd = s2chrs[sp]; //species data object
	if(!spd){
	    spd = {
		species : sp,
		chrs : [],
		maxLen : -1
	        };
            // make sure the aName species goes into pos 0 and bName species into 1.
            species[ (sp === aName) ? 0 : 1 ] = spd;
	    s2chrs[sp] = spd;
	}
	dd = {
	    id : chr,
	    length : len,
	    scale : d3.scale.linear().domain([1,len]).range([0,height]),
	    species : sp
	    };
	spd.chrs.push(dd);
	spd.maxLen = Math.max(spd.maxLen, len);
	maxLen = Math.max(maxLen, spd.maxLen); 
    });
    yscale = d3.scale.linear().domain([1,maxLen]).range([0,height]);
    species.forEach(function(spd){
	spd.xscale = d3.scale.ordinal()
	  .domain(spd.chrs.map(function(x){return x.id;}))
	  .rangePoints([0,width], 2);
	spd.yscale = d3.scale.linear()
	  .domain([1,spd.maxLen])
	  .range([0,height]);
	spd.cscale = d3.scale
            .category20().domain(spd.chrs.map(function(x){return x.id;}));
	spd.chrs.forEach(function(c){
	    var sc = d3.scale.linear().domain([1,c.length]).range([0, yscale(c.length)]);
            c.brush = d3.svg.brush().y(sc)
               .on("brushstart",brushstart)
               .on("brushend",brushend);
	  });
    });

    colors = species[0].cscale;

    d3.selectAll('a.sblock').remove();

    redraw();
}

// ---------------------------------------------

function setZoomLabel(scale){
  if(!arguments.length)
      scale = d3.event.scale;
  var ys = yscale.range();
  var yr = (ys[1]-ys[0])/scale;
  yViewSize = Math.round(yr * 1000000)
  text.text(formatLength(yViewSize));
}

function zoom() {
  svg.attr("transform", "translate(" + d3.event.translate + ")scale(" + d3.event.scale + ")");
  magnification = d3.event.scale;
  setZoomLabel();
  drawBlocks();
}

function resetView(){
  zoomer.scale(1).translate([0,0]);
  svg.attr("transform", "translate(0,0)scale(1)");
  magnification = 1;
  setZoomLabel(1);
  drawBlocks();
}

function setTitle(bspecies, cspecies) {
    // Species label
    //
    var b;
    var slabel = svg.selectAll("#specieslabel").data([bspecies]);
    slabel.enter().append('text')
        .attr("id", "specieslabel")
        .attr("x", 0.3* width)
	.attr("y", -20)
	;
    slabel.text(["Genome = " + bspecies, " Color = " + cspecies]);

    d3.select('#bFlipGenome').text("Genome: "+bspecies);
    d3.select('#bFlipColors').text("Colors: "+cspecies);

}


function blockClicked () {
    let alt = d3.event.altKey;
    let b = d3.select(this);
    let ids = b.data()[0].ids;
    //
    let formName = alt ? 'mmLinkForm' : 'mgiLinkForm';
    let inpName  = alt ? 'externalids' : 'ids';
    let joinChar = alt ? ',' : ' ';
    let form = d3.select(`form[name="${formName}"]`);
    let input = form.select(`input[name="${inpName}"]`);
    input[0][0].value = ids.join(joinChar);
    //
    form[0][0].enctype = alt ? '' : 'multipart/form-data';
    form[0][0].submit()
}

function redraw () {
    drawChrs( species[spIndex], species[cIndex], colors);
}

function drawBlocks(){
    let s= 'a';
    d3.selectAll('a.sblock rect')
        .attr("height", function(d){
             return Math.max(minRectHeight/magnification, yscale(d[s+'End']-d[s+'Start']+1));
        })
}

function drawChrs(spd, spd2, colors){
    var ccells;
    var slabel = null;
    //var labels;
    var brushes;
    var blks;
    var xf = function(d){return bwidth+spd.xscale(d.id);};
    var s,s2;

    setTitle(spd.species, spd2.species);

    // Chromosome labels
    //
    labels = svg.selectAll('.chrlabel')
      .data(spd.chrs, function(x){return x.id;});
    labels.exit().transition().duration(dur)
      .attr('y', height)
      .remove();
    labels.enter().append('text')
      .attr('class','chrlabel')
      .attr('x', xf)
      .attr('y', height)
      .on('click', chrClicked);
    labels
      .text(function(d){return d.id;})
      .transition().duration(dur)
      .attr('x', xf)
      .attr('y', -2) ;

    // Chromosome backbones (lines)
    //
    ccels = svg.selectAll('line.chr')
      .data(spd.chrs, function(x){return x.id;});
    ccels.exit().transition().duration(dur)
      .attr("y1", height)
      .attr("y2", height)
      .remove();
    ccels.enter().append('line')
      .attr('class','chr')
      .attr("x1", xf)
      .attr("y1", height)
      .attr("x2", xf)
      .attr("y2", height)
      ;
    ccels.transition().duration(dur)
      .attr("x1", xf)
      .attr("y1", 0)
      .attr("x2", xf)
      .attr("y2", function(d){return yscale(d.length);})
	;

    // Synteny blocks -----------------------------
    //
    // Block data
    s  = ["a","b"][spIndex];
    s2 = ["a","b"][cIndex];
    let blocksToDraw = allBlocks;
    if(facets.length > 0){
        blocksToDraw = blocksToDraw.filter(function(b) {
            let vals = facets.map( f => facetFuncs[f](b) );
            let v = vals.reduce( (acc,cv) => acc && cv, true );
            return v;
        });
    }
    blks = svg.selectAll('a.sblock')
      .data(blocksToDraw, function(x){return x.name;});


    // Block DOM structure
    blks.enter()
       .append('a')
       .attr('class','sblock')
       .attr('target','_blank')
         .append('rect')
	 .attr('class', 'sblock')
         .on("click", blockClicked)
           .append('title')    
	   ;

    //
    blks.exit().remove();

    // Block tooltips
    blks.select('title').text(function(d,i){
         return `A=${fmtLoc(d.aChr,d.aStart,d.aEnd)}(${formatLength(d.aLength)})`
            +   ` B=${fmtLoc(d.bChr,d.bStart,d.bEnd)}(${formatLength(d.bLength)})`
            +   `\n${d.ids.length} gene${d.ids.length > 1 ? 's':''}`
            +   `\n${d.ids.slice(0,5).join(' ')}`
            +   (d.ids.length > 5 ? "..." : "");
     })

    // Transition effects
    blks.select('rect')
     .classed("inverted", x => x.ori == "-")
     .attr("fill", function(d){ return colors(d[s2+'Chr']); })
     .transition().duration(dur)
     .attr("x", function(d){
         let x = spd.xscale(d[s+'Chr']) + (d.ori=="+"?bwidth:0);
         if (isNaN(x)) throw "Bad coordinate.";
         return x;
     })
     .attr("y", function(d){return yscale(d[s+'Start']);})
     .attr("width",bwidth) 
     .attr("height", function(d){
         return Math.max(minRectHeight/magnification, yscale(d[s+'End']-d[s+'Start']+1));
     })
     ;

     // Map functions
     blks.each(function(d){
       var s = spd.species === aName ? "a" : "b";
       var s2= spd2.species === aName ? "a" : "b";
       d.map.domain([d[s+'Start'],d[s+'End']]);
       d.map.range(d.ori=='+'?[d[s2+'Start'],d[s2+'End']]:[d[s2+'End'],d[s2+'Start']]);
     });

    // Brushes
    brushes = svg.selectAll("g.brush").data(spd.chrs, function(x){return x.id;});
    brushes.exit().remove();
    brushes.enter().append('g').attr('class','brush')
        .each(function(d){d3.select(this).call(d.brush);})
	.selectAll('rect')
	 .attr('width',10)
	 .attr('x', cwidth/4)
	 ;
    brushes
        .attr('transform', function(d){return 'translate('+spd.xscale(d.id)+')';})
        .each(function(d){d3.select(this).call(d.brush);})
	;
}

function brushstart(c){
    clearBrushes(c.brush);
    d3.event.sourceEvent.stopPropagation();
}

function brushend(c){
    if(c.brush.empty()) return;
    var xtnt = c.brush.extent();
    console.log("Brush end.", c.id, xtnt)
}

function clearBrushes(except){
    d3.selectAll('.brush').each(function(c){
	if(c.brush !== except){
	    c.brush.clear();
	    c.brush(d3.select(this));
	}
    });
}

function fmtLoc(c,s,e){
    s = Math.floor(s);
    e = Math.floor(e);
    return (s<=e ? c+":"+s+".."+e : c+":"+e+".."+s);
}

function clearSelections(){
    clearBrushes();
    aSelected = [];
    bSelected = [];
    fadeEm(aSelected, bSelected);
}

function flipGenome(){
    clearBrushes();
    spIndex = 1-spIndex;
    redraw();
}

function flipColors(){
    cIndex = 1-cIndex;
    var s2 = ["aChr","bChr"][cIndex];
    d3.selectAll('.sblock')
      .transition().duration(1000)
      .attr('fill', function(d){ return colors(d[s2]); });
    setTitle( species[spIndex].species, species[cIndex].species );
}

function chrClicked(d){
    let lst = d.species === aName ? aSelected : bSelected;
    let i = lst.indexOf(d.id);
    if(i==-1)
        lst.push(d.id);
    else
        lst.splice(i,1);
    fadeEm(aSelected, bSelected);
}

function fadeEm(aChrs, bChrs){
    svg.selectAll('.sblock').classed('fade', function(d){
	return (aChrs.length>0 && aChrs.indexOf(d.aChr) == -1) || (bChrs.length>0 && bChrs.indexOf(d.bChr) == -1);
        });
}

