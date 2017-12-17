var mmServer = "www.mousemine.org";

var mousemine   = new intermine.Service({root: mmServer+'/mousemine/service'});
var margin = {top: 35, right: 10, bottom: 20, left: 30};
var outerWidth = 800;
var outerHeight = 400;
var width = outerWidth - margin.left - margin.right;
var height = outerHeight - margin.top - margin.bottom;
var zoomer = d3.behavior.zoom().on("zoom", zoom);
var labels = null;
var maxBlockSize = -1; // Max size in bp of synteny blocks to draw. -1 for unlimited
var minRectHeight = 0; // Min size in px for drawing a block. 0 for no minimum
var magnification = 1.0;
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

d3tsv("./output/strainList.tsv").then(function(strains) {
    let s1 = initOptList("#aGenome", strains, s=>s.strain)
    let s2 = initOptList("#bGenome", strains, s=>s.strain)
    d3.select("#go")
        .on("click", go);
});

function go () {
    let s1 = d3.select("#aGenome")[0][0].value;
    let s2 = d3.select("#bGenome")[0][0].value;
    console.log("GO!", s1, s2);
    let makefn = (a,b) => `./output/${a}-${b}.tsv`;
    Promise.all([
        d3tsv(makefn(s1,s1)),
        d3tsv(makefn(s2,s2)),
        d3tsv(makefn(s1,s2))
    ]).then(function(data) {
        let aChrs = data[0].map(c => [ c.aChr, c.aLength, s1 ]);
        let bChrs = data[1].map(c => [ c.bChr, c.bLength, s2 ]);
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
            mchr   : k.aChr,
            mstart : k.aStart,
            mend   : k.aEnd,
            hchr   : k.bChr,
            hstart : k.bStart,
            hend   : k.bEnd,
            map    : d3.scale.linear().clamp(true)
          });
        });
        window.sblocks = bks;
        processChrs( allChrs );
    });
}

function setMaxBlockSize(n) {
    maxBlockSize = n;
    redraw();
}

function setMinRectHeight(n) {
    minRectHeight = n;
    redraw();
}

// ---------------------------------------------

function setText(scale){
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
  setText();
}

function resetView(){
  zoomer.scale(1).translate([0,0]);
  svg.attr("transform", "translate(0,0)scale(1)");
  setText(1);
}

function processChrs(data) { 
    // data = Data for each chromosome in order by id
    // (chrs from different species are interleaved)
    //
    // Data = list of: primaryIdentifier, length, organism.shortName
    //
    var sp, chr, len;
    var species = [];
    var s2chrs = {};
    var spd, dd; // spd=species data, dd=data for 1 chr
    var sx,sy;
    var data2 = [];
    var maxLen = 0; // for any chromosome of any species
    data.forEach(function(d){
	sp = d[2];
	chr = d[0];
	len = +d[1];
	if(len < 100000) return; // skip contigs
	spd = s2chrs[sp];
	if(!spd){
	    spd = {
		species : sp,
		chrs : [],
		maxLen : -1
	        };
	    species.push(spd);
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
    window.yscale = d3.scale.linear().domain([1,maxLen]).range([0,height]);

    window.species = []; // list of species data
    window.spIndex = 0; // index of species for backbone
    window.cIndex = 0; // index of species for coloring
    window.colors=null; // color map

    species.forEach(function(spd){
	spd.xscale = d3.scale.ordinal()
	  .domain(spd.chrs.map(function(x){return x.id;}))
	  .rangePoints([0,width], 2);
	spd.yscale = d3.scale.linear()
	  .domain([1,spd.maxLen])
	  .range([0,height]);
	spd.cscale = d3.scale.category20c().domain(spd.chrs.map(function(x){return x.id;}));
	spd.chrs.forEach(function(c){
	    var sc = d3.scale.linear().domain([1,c.length]).range([0, window.yscale(c.length)]);
	    c.brush = d3.svg.brush().y(sc)
	      .on("brushstart",brushstart)
	      .on("brushend",brushend)
	      ;
	  });
	window.species.push(spd);
    });

    window.colors = window.species[0].cscale;

    redraw();
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

function redraw () {
    drawChrs( species[spIndex], species[cIndex], window.colors);
}

function drawChrs(spd, spd2, colors){
    var ccells;
    var brushes;
    var slabel = null;
    //var labels;
    var blks;
    var cwidth = 20;
    var dur = 1500;
    var bwidth = cwidth/2;
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
      .attr("y2", function(d){return window.yscale(d.length);})
	;

    // Synteny blocks
    //
    s = ["m", "h"][spIndex];
    s2 = ["m","h"][cIndex];
    let blocksToDraw =
        maxBlockSize === -1 ? window.sblocks
            : window.sblocks.filter(b => b.hend - b.hstart+1 <= maxBlockSize);
    blks = svg.selectAll('a.sblock')
      .data(blocksToDraw, function(x){return x.name;});
    blks.enter()
       .append('a')
       .attr('class','sblock')
       .attr('target','_blank')
         .append('rect')
	 .attr('class', 'sblock')
           .append('title')    
	   ;

    blks.exit().remove();

    blks.select('title').text(function(d,i){
         return d.ori + ' A'+fmtLoc(d.mchr,d.mstart,d.mend)+' B'+fmtLoc(d.hchr,d.hstart,d.hend);
     })
    blks.select('rect')
     .classed("inverted", x => x.ori == "-")
     .transition().duration(dur)
     .attr("x", function(d){return spd.xscale(d[s+'chr']) + (d.ori=="+"?bwidth:0);})
     .attr("y", function(d){return window.yscale(d[s+'start']);})
     .attr("width",bwidth) 
     .attr("height", function(d){
         return Math.max(minRectHeight/magnification, window.yscale(d[s+'end']-d[s+'start']+1));
     })
     .attr("fill", function(d){ return colors(d[s2+'chr']); })
     ;
    blks.each(function(d){
     var s = spd.species[0].toLowerCase();
     var s2= spd2.species[0].toLowerCase();
     d.map.domain([d[s+'start'],d[s+'end']]);
     d.map.range(d.ori=='+'?[d[s2+'start'],d[s2+'end']]:[d[s2+'end'],d[s2+'start']]);
     });

    /*
    // Brushes
    //
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
        */

}

function fmtLoc(c,s,e){
    s = Math.floor(s);
    e = Math.floor(e);
    return (s<=e ? c+":"+s+".."+e : c+":"+e+".."+s);
}

function brushstart(c){
    clearBrushes(c.brush);
    d3.event.sourceEvent.stopPropagation();
}

function brushend(c){
    if(c.brush.empty()) return;
    var s = c.species[0].toLowerCase();
    var s2 = (s=='m'?'h':'m');
    var s2n = s2=='m'?'M. musculus' : 'H. sapiens';
    var xtnt = c.brush.extent();
    var tx1 = (s=="m"?10090:9606);
    var tx2 = (s=="m"?9606:10090);
    var r;
    var mmOverlapQuery = {
      from: 'SequenceFeature',
      select: [
	'organism.shortName',
	'primaryIdentifier',
	'symbol',
	'chromosomeLocation.locatedOn.primaryIdentifier',
	'chromosomeLocation.start',
	'chromosomeLocation.end',
	'chromosomeLocation.strand'
      ],
      where : [],
      sortOrder: [["organism.shortName", "ASC"], ["symbol","ASC"]]
    };

    // A
    mmOverlapQuery.where.push(["primaryIdentifier", "is not null"]);

    // (B and C)
    mmOverlapQuery.where.push(["organism.taxonId", "=", tx1]);
    r = fmtLoc(c.id, Math.floor(xtnt[0]), Math.floor(xtnt[1]));
    mmOverlapQuery.where.push(["chromosomeLocation","OVERLAPS",[r]]);


    var rs = [];
    sblocks.forEach(function(b){
	if(c.id == b[s+'chr'] && xtnt[0]<=b[s+'end'] && xtnt[1] >= b[s+'start']){
	    rs.push( fmtLoc(b[s2+'chr'], Math.floor(b.map(xtnt[0])), Math.floor(b.map(xtnt[1]))) );
	}
    });
    if(rs.length > 0){
	// (D and E)
	mmOverlapQuery.constraintLogic = "A and ((B and C) or (D and E))";
	mmOverlapQuery.where.push(["organism.taxonId", "=", tx2]);
	mmOverlapQuery.where.push(["chromosomeLocation", "OVERLAPS",rs]);
    }

    if(mmOverlapQuery){
	console.log(mmOverlapQuery);
	mousemine.rows(mmOverlapQuery).then(function(rows){ 
	  var rs = d3.select('#genestbl').selectAll('tr').data(rows);
	  rs.enter().append('tr');
	  rs.exit().remove();
	  var tds = rs.selectAll('td').data(function(d){return d;});
	  tds.enter().append('td');
	  tds.text(function(d){return d;});
	});

	window.currOverlapQuery = mmOverlapQuery;
	d3.select('#bToMouseMine').attr('disabled',null);
    }

}

function runOverlapQuery(){
    var mmOverlapQuery = window.currOverlapQuery;
    if(! mmOverlapQuery ) return;

    /* Run the query and get back rows. Dump symbols to page. */
    console.log(mmOverlapQuery);
    // Turn query into a link to execute at MouseMine. Open in new page.
    var q = new intermine.Query(mmOverlapQuery).toXML();
    q=q.replace('query','query model="genomic" ');
    q = mmRunQueryLink.replace('@@@@', encodeURIComponent(q));
    window.open(q);
}

function clearBrushes(except){
    window.currOverlapQuery = null;
    d3.select('#bToMouseMine').attr('disabled',true);

    d3.selectAll('.brush').each(function(c){
	if(c.brush !== except){
	    c.brush.clear();
	    c.brush(d3.select(this));
	}
    });
}

function clearSelections(){
    clearBrushes();
    window.mselected = [];
    window.hselected = [];
    fadeEm(window.mselected, window.hselected);
}

function flipGenome(){
    clearBrushes();
    spIndex = 1-spIndex;
    redraw();
}

function flipColors(){
    cIndex = 1-cIndex;
    var s2 = ["m","h"][cIndex];
    d3.selectAll('.sblock')
      .transition().duration(1000)
      .attr('fill', function(d){ return colors(d[s2 + 'chr']); });
    setTitle( species[spIndex].species, species[cIndex].species );
}

window.mselected=[]; 
window.hselected=[];

function chrClicked(d){
    var s = d.species[0].toLowerCase() + 'selected';
    var i = window[s].indexOf(d.id);
    if(i==-1)
        window[s].push(d.id);
    else
        window[s].splice(i,1);
    fadeEm(window.mselected, window.hselected);
}

function fadeEm(mchrs, hchrs){
    svg.selectAll('.sblock').classed('fade', function(d){
	return (mchrs.length>0 && mchrs.indexOf(d.mchr) == -1) || (hchrs.length>0 && hchrs.indexOf(d.hchr) == -1);
        });
}

