var wTotal=1010, wTreeSvg=432, wPheno=140, wSvg=wTotal-wTreeSvg-wPheno, spTip=wTreeSvg+wPheno+55;
var hUnit=14, hBlock=5, hTick=4, lTri=8, hSvg=45, base_tick=17, base=base_tick+14;
var wCell=40, wTbl, hPos=57, snpPos, snpSeq, snpNum, vSnpLine=9, eSnpLine=14;

var contig = {'NC_003197':'Chromosome','NC_003277':'Plasmid'}, isolates, ids=['ref'], brch2pos, pos2brch;
var bid_N2O = {2:3,3:17,4:2,5:5,6:16,7:7,8:9,9:10,10:8,11:12,12:15,13:14,14:13,15:11,16:6,17:4},
	bid_O2N = {2:4,3:2,4:17,5:5,6:16,7:7,8:10,9:8,10:9,11:15,12:11,13:14,14:13,15:12,16:6,17:3};

var r=25, genomeGap=10*r;

$(document).ready(function(){
    $('#phenoTbl').css("width", wPheno+'px');
    $('#showSnp, #svgOrf').css("width", wSvg+'px');

    $.ajax({
        async: false,
        url: "js-css/data.json",
        dataType: "json",
        success: function(data) {
            isolates=data.iso;
            snpPos=data.snp; snpSeq=data.seq; snpNum=data.snpNum;
            brch2pos=data.b2p; pos2brch=data.p2b;
            drawORF(data.length, data.orf)
        }
    });

    $.ajax({
        async: false,
        url: "js-css/tree.json",
        dataType: "json",
        success: function(data) {
        	$.each(data.ids, function(n,id){ ids.push(id)});
        	draw_tree('#trees', data) }
    });

    phenoTbl();

	wTbl = wCell * snpPos.length;
	snpTbl();

    $('#showSnp').height(hPos + hUnit*ids.length + 24);

    $('#snpTbl').css("width", wTbl+'px');
    $('#snpTbl td').css("width", wCell);
    $('#snpTbl td, #phenoTbl td').css("line-height", (hUnit-1)+'px')

    $('#closeA').on("click", function(){ $('#answer').hide(); return false});
    $('#closeM').on("click", function(){ $('#method').hide(); return false});

//	draw_legORF('#legORF');
});

function showA(){ $('#answer').show() }
function showM(){ $('#method').show() }


function initMap() {
	var city = {'1_1':[37.540725,-77.436048],
				   '2_1':[40.712784,-74.005941],
				   '2_2':[42.652579,-73.756232],
				   '2_3':[42.443961,-76.501881],
			   	   '2_4':[41.035935,-71.954515]};

	var map = new google.maps.Map(document.getElementById('map'), {
		disableDefaultUI:true,
		scaleControl:true,
		center: new google.maps.LatLng(city['2_1'][0],city['2_1'][1]),
		zoom:5
    });

	$.each(city, function(ss,obj){
		var arr = ss.split('_'),
			col = arr[0]==1? "rgba(0,0,255,0.4)" : "rgba(255,0,0,0.4)";
	    var marker = new google.maps.Circle({
					center: new google.maps.LatLng(obj[0],obj[1]),
					radius: 30000,
					strokeColor: col,
					strokeWeight: 1,
					fillColor: col,
					fillOpacity: 1,
					map: map,
					title: 'state ' + arr[0] + ', site ' + arr[1]
  		});

		google.maps.event.addListener(marker,'mouseover',function(){
             this.getMap().getDiv().setAttribute('title',this.get('title'))});
		google.maps.event.addListener(marker,'mouseout',function(){
             this.getMap().getDiv().removeAttribute('title')});
		google.maps.event.addListener(marker,'click',function(){
			$('#phenoTbl td').css('background-color', '');
			$('.s_'+ss).css('background-color', col) })
	})
}

var p0;
function snpTbl() {
    var svg = d3.select('#position').append("svg").attr("width", wTbl).attr("height", hPos),
    	svga = svg.append("svg:a").attr("xlink:href", "#");
    $.each(snpPos, function(i,obj){
		svga.append("text")
			.attr("id", 'p_' + (i+1))
			.attr("class",obj[1]? (obj[2]? 'posS' : 'posN') : 'pos')
		    .attr("x", 0).attr("y", 0)
	    	.attr("transform", "translate(" + (wCell*i+23-(obj[1]? 5*(3-obj[1]-(obj[2]? 0 : 1)) : 0)) + "," + (hPos-1) + "),rotate(-90)")
	    	.text(numberWithCommas(obj[0]))
			.on("click", function(){ highLight(i+1); scrollWin(i+1) })
    });

    var table='<table id=' + 'snpTbl' + '>';
    $.each(ids, function(n,gId){
 		table += '<tr>';
 		var pMin = 0;
 		$.each(contig, function(cont,name){
 			var thisSeq = snpSeq[cont][gId],
 				posN = snpNum[cont];
			$.each(thisSeq, function(i,s){
				var p=snpPos[pMin+i], c='';
				if (gId=='ref' || !p[1]){ c = s[0] }
				else { for (var j=1; j<=3; j++){ c += p[1]==j? s[0] : '&nbsp' } }
				table += '<td>' + c;

				if (p[2]){
					var a = p[2];
					table += ' <span class=Lewis_' + a[s[1]?1:0] + '>' + (gId=='ref'? a[0] : (s[1]? a[1]:'.')) + '</span>'
				}
 			});
			table += '</td>';
			pMin += posN
		});
		table += '</tr>'
	});

    table+='</table>';
    $("#showSnp").append(table)
}


var lowerYbound=11;
function draw_tree(container, data) {
    var nodes = data.nodes,
    	w_tree = data.w_tree,
    	unit,
    	svg_h = hUnit*(ids.length-1)+4;

    var wTreeTip=57,
    	wTree = wTreeSvg-wTreeTip,
    	lowerXbound=1,
    	gap=5,
    	r = wTree/w_tree;

    //unit
    $.each(nodes, function(index,node){
    	if (node.id=='SRR2352185'){ unit = node.branch_length*r/brch2pos[node.id].length; return }
    });

    var svg = d3.select(container).append("svg").attr("width", wTreeSvg).attr("height", svg_h);

    var line_tree = svg.append("g").style("stroke", "black"),
    	line_a = line_tree.append("svg:a").attr("xlink:href", "#"),
    	lineDotted = svg.append("g").attr("class", "lineDotted"),
    	text_a = svg.append("svg:a").attr("xlink:href", "#");

    // draw tree:
	var brNum=0;
    $.each(nodes, function(index,node){
		var xx = node.xcoord*r + lowerXbound,
			yy = node.ycoord*hUnit + lowerYbound,
			bl = node.branch_length*r;

	    if (!node.is_Leaf) { brNum++ }

	    if (bl){
	    	var id = node.is_Leaf? node.id : brNum;
			line_tree.append("line")
				  .attr("x1", xx).attr("x2", xx-bl).attr("y1", yy).attr("y2", yy)
				  .attr("id", 'br_'+id);

			var pos = brch2pos[node.is_Leaf? node.id : bid_N2O[brNum]],
				xstart = xx - bl + unit/2;

			$.each(pos, function(n, p){
				var obj = snpPos[p-1];
				if (!obj) { return }
				var x = xstart + unit*n;
			    line_a.append("line")
	    			  .attr("id", 'm_'+p)
	    			  .attr("x1", x).attr("x2", x)
	    			  .attr("y1",yy-vSnpLine/2).attr("y2", yy+vSnpLine/2)
		    		  .attr("class", obj[1]? (obj[2]? "snpLineS" : "snpLineN") : "snpLine")
		    		  .on("click", function(){ highLight(p); scrollWin(p,1) })
		    		  .on("mouseover", function() {
					    	d3.select("#tipBrch")
					    	  .style("left", (x-25)+"px")
					    	  .style("top", yy-8 + "px")
					    	  .html(numberWithCommas(obj[0])).classed("hidden", false) })
					  .on("mouseout", function() { d3.select("#tipBrch").classed("hidden", true) })
			})
		}

//		svg.append("circle").attr("cx", xx).attr("cy",yy).attr("r", 1.5)

		if (node.is_Leaf) {
	    	text_a.append("text")
	    		  .attr("x",wTree+lowerXbound+gap).attr("y",yy+3)
	    		  .text(isolates[node.id].id)
	    		  .on("click", function(){ window.open("http://trace.ncbi.nlm.nih.gov/Traces/sra/?run="+node.id) })
	    	lineDotted.append("line")
					  .attr("x1", xx+4).attr("x2", wTree+lowerXbound)
					  .attr("y1", yy).attr("y2", yy)
		} else {
	    	line_tree.append("line")
	    			 .attr("x1", xx).attr("x2", xx)
	    			 .attr("y1", node.descendent_ycoord[0]*hUnit+lowerYbound)
	    			 .attr("y2", node.descendent_ycoord[1]*hUnit+lowerYbound);
	    	if (brNum==1) {return}
/*	    	svg.append("circle").attr("cx", xx).attr("cy",yy).attr("r", 2)
	    		   .attr("id", 'node_'+brNum)
	    		   .on("mouseover", function() {
				    	d3.select("#tipBrch")
				    	  .style("left", (xx-25)+"px")
				    	  .style("top", yy+10 + "px")
				    	  .html(this.id).classed("hidden", false) })
					.on("mouseout", function() { d3.select("#tipBrch").classed("hidden", true) });*/
		}
    });

    //draw scale:
    var scaleY = lowerYbound + 334, scaleX = 16;
    line_tree.append("line").attr("x1", scaleX).attr("x2", scaleX+unit).attr("y1", scaleY).attr("y2", scaleY);
    line_tree.append("line").attr("x1", scaleX).attr("x2", scaleX).attr("y1", scaleY-1.5).attr("y2", scaleY+1.5);
    line_tree.append("line").attr("x1", scaleX+unit).attr("x2", scaleX+unit).attr("y1", scaleY-1.5).attr("y2", scaleY+1.5);

    svg.append("text").attr("id",'treeLg').attr("x", scaleX+unit+7).attr("y", scaleY+3).text('1 substitution')
}


function phenoTbl() {
	var trs=[];
	$.each(ids, function(n,acc) {
		var tr = '<tr>';
		if (acc=='ref'){
			tr += '<td>&nbsp</td><td></td><td>'
		} else {
			var obj = isolates[acc]
			tr += (obj.source=='outbreak' || obj.source=='human'? '<td class='+obj.source+'>' : '<td>') + obj.source + '</td>';
			if (obj.site) {
				var cls = 's_' + obj.state + '_' + obj.site;
				tr += '<td class=' + cls + '>' + obj.state + '</td><td class=' + cls + '>' + obj.site
			} else {
				tr += '<td>' + (obj.state? obj.state : '') + '</td><td>'
			}
		}
		tr += '</td></tr>';
		trs.push(tr)
	});
	$('#phenoTbl').append(trs.join('')).append('<tr><th>source</th><th>state</th><th>site</th></tr>').show()
}


function drawORF(orfLen, orfs) {
	var totalL=0;
	$.each(contig, function(cont,name){
		if (totalL){ totalL += genomeGap }
		totalL += orfLen[cont]
	});

    var svg = d3.select('#svgOrf') .append("svg") .attr("width", totalL/r+1) .attr("height", hSvg);
    var axis = svg.append("g").attr("id", 'axis'),
		line = svg.append("g").attr("class", "baseline"),
	    svga = svg.append("svg:a").attr("xlink:href", "#");

	var xscale = d3.scale.linear().domain([0, 100]).range([0, 100/r]),
		x0=0,
		pMin=0;

	$.each(contig, function(cont,name){
		var lGenome = orfLen[cont],
			posN = snpNum[cont];

		//draw axis
    	drawAxis(axis, 200, 10, 1000, r, x0, lGenome, name);

		//draw base line
		line.append("line").attr("x1",xscale(x0)).attr("x2",xscale(x0+lGenome)).attr("y1",base).attr("y2",base);
		line.append("line").attr("x1",xscale(x0)).attr("x2",xscale(x0)).attr("y1",base-3).attr("y2",base+3);
		line.append("line").attr("x1",xscale(x0+lGenome)).attr("x2",xscale(x0+lGenome)).attr("y1",base-3).attr("y2",base+3);

		//draw orfs
		$.each(orfs[cont], function(n, obj){
	    	var end = xscale(x0 + obj.end*1),
	   			start = xscale(x0 + obj.end*1 - obj.length);

		   	var thisOrf;
		   	if (obj.rna){
				thisOrf = svga.append("rect")
				   		  .attr("x", obj.length>0? start : end)
			   			  .attr("y", base-hBlock/2+0.5)
			   			  .attr("width", Math.abs(end-start))
				   		  .attr("height", hBlock-1)
				   		  .attr("class", "rna")
	   		} else {
		    	var xv = [start, end-lTri*(obj.length>0? 1 : -1), end, end-lTri*(obj.length>0? 1 : -1), start],
		    		yv = [base-hBlock/2, base-hBlock/2, base, base+hBlock/2, base+hBlock/2];
				var points = [];
				for (var i=0; i<xv.length; i++) { points.push(xv[i] + ',' + yv[i]) }

			   thisOrf = svga.append("polygon")
						 .attr("points", points.join(' '))
						 .attr("id", obj.locus)
						 .attr('end', obj.end)
						 .attr('len', obj.length)
						 .attr("class", obj.pseudo ? "pseu" : "orf")
	   		}
		   	var tipText = obj.locus;
		   	if (obj.sym) { tipText += '&nbsp;&nbsp;<em>' + obj.sym + '</em>' }
	   		thisOrf.on("mouseover", function(){
						var xMouse = d3.event.pageX + $("#svgOrf").scrollLeft() - wTreeSvg - wPheno-40;
						d3.select("#tooltip").style("left", xMouse+"px")
											 .style("top", (base-21) + "px")
											 .classed("hidden", false)
											 .html(tipText)
						})
					.on("mouseout", function(){ d3.select("#tooltip").classed("hidden", true) })
					.on("click", function(){ window.open("http://www.ncbi.nlm.nih.gov/gene/?term="+obj.locus) })
		});

		//draw position lines
		for (n=pMin; n<pMin+posN; n++){
			var obj = snpPos[n];
		    svg .append("line")
	    		.attr("id", 'l_'+(n+1))
	    		.attr("x1", xscale(x0+obj[0])).attr("x2", xscale(x0+obj[0]))
		    	.attr("y1",base-vSnpLine).attr("y2", base+vSnpLine)
		    	.attr("class", obj[1]? (obj[2]? "snpLineS" : "snpLineN") : "snpLine")
		}
		x0 += lGenome + genomeGap;
		pMin += posN
	})
}

function drawAxis(fig, unit, sub, k, r, x0, lGenome, ct) {
    fig.append("line").attr("x1", (x0+1)/r).attr("x2", (x0+lGenome)/r).attr("y1", base_tick).attr("y2", base_tick);
    fig.append("line").attr("x1", (x0+1)/r).attr("x2", (x0+1)/r).attr("y1", base_tick-hTick/2).attr("y2", base_tick+hTick).style("stroke-width", "1.6px");
//    fig.append("text").attr("x",(x0+1)/r).attr("y", base_tick-5).text(1).style("text-anchor", "start");
    fig.append("text").attr("x",(x0+1)/r).attr("y", base_tick-9).text(ct).style({"text-anchor":"start", "fill":"black"});
    for (i=1; i<=parseInt(lGenome/unit); i++){
		var x = (x0+i*unit)/r;
		fig.append("line").attr("x1",x).attr("x2",x).attr("y1",i%5? base_tick : base_tick-hTick/2).attr("y2",base_tick+(i%5? hTick/2 : hTick));
		if (i%sub==0) { fig.append("text").attr("x", x).attr("y", base_tick-5).text(numberWithCommas(i*unit/k) + (k==1000? 'k' : (k==1000000? 'M' : ''))) }
    }
}

function numberWithCommas(x) { return x.toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",") }

function highLight(p) {
	d3.selectAll('#position text').classed("bigger", false);
	d3.selectAll('#trees line').classed("changedLine", false);
	if (p0) {
		d3.selectAll('#m_'+p0).classed("lineWidest",false);
		d3.select('#l_'+p0).attr("y1", base-vSnpLine).attr("y2", base+vSnpLine).classed("lineWider",false);
		$('#snpTbl td:nth-child('+p0+')').css('background-color', '');
	}
	p0 = p;

	d3.select('#p_'+p).classed("bigger", true);
	d3.selectAll('#m_'+p).classed("lineWidest",true);
	d3.select('#l_'+p).attr("y1", base-eSnpLine).attr("y2", base+eSnpLine).classed("lineWider",true);
	$('#snpTbl td:nth-child('+p+')').css('background-color', '#FFC');
	$.each(pos2brch[p], function(n,b){
		d3.select('#br_' + ($.isNumeric(b)? bid_O2N[b] : b)).classed("changedLine",true)
	})
}

function scrollWin(p, clickOrf) {
	var t1 = d3.select('#l_'+p).attr("x1"),
		t2 = d3.transform(d3.select('#p_'+p).attr("transform")).translate[0];

	if (clickOrf){ $("#showSnp").scrollLeft(t2 - wSvg/2) }

	var pText = t2-$('#showSnp').scrollLeft();
	$("#svgOrf").scrollLeft(t1 - (pText-4))
}


/*
function draw_legORF(container) {
    var svg_w = 385,
    	svg_h = 19,
    	base = 8,
		xPos = [5,98,171,240,313],
		recCls = ['oidLeg', 'oidLeg', 'oidLeg', 'oidSingle', 'oidSingle'],
		lTxt = ['core ortholog', 'homolog', 'paralog', 'singleton', 'pseudogene'],
		w = hBlock*3+1;

    var svg = d3.select(container).append("svg").attr("width", svg_w).attr("height", svg_h);
    var txt_leg = svg.append("g").attr("class", "txt_leg");

    $.each(xPos, function(index,x) {
		var rect = svg.append("rect")
			.attr("x", x).attr("y", base)
			.attr("width", w)
			.attr("height", hBlock)
			.attr("class", recCls[index]);
 		if (index>0) {
 			rect.style("fill", "none");
 			if (index==2 || index==4) { rect.style({"stroke-dasharray":"2,1", "stroke-width":"0.5px"}) }
 		}

 		txt_leg.append("text").attr("x", x+w+4).attr("y", base+5).text(lTxt[index])
    })
}
*/