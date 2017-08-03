htmlTitle <- function(outputHtmlFile,hostName,organism,geneSetNet,geneSetDAG){
	 cat("<html>\n",file=outputHtmlFile)
	 cat("<head>\n",file=outputHtmlFile,append=TRUE)
	 cat("<title>WebGestalt (WEB-based GEne SeT AnaLysis Toolkit)</title>\n",file=outputHtmlFile,append=TRUE)
	 cat('<script type="text/javascript" src="',file.path(hostName,"js","jquery.js"),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
	 cat('<script type="text/javascript" src="',file.path(hostName,"js","jquery-ui.js"),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
	 cat("<script>(function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
  (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
  m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
  })(window,document,'script','//www.google-analytics.com/analytics.js','ga');
  ga('create', 'UA-42098540-1', 'auto');
  ga('send', 'pageview');
</script>\n",file=outputHtmlFile,append=TRUE,sep="")
	 if(organism!="others"){
			 cat('<script type="text/javascript" src="',file.path(hostName,"js","changeColor.js"),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
		   cat('<script type="text/javascript" src="',file.path(hostName,"js","tableExport.js"),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
		   cat('<script type="text/javascript" src="',file.path(hostName,"js","jquery.base64.js"),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
		   cat('<script type="text/javascript" src="',file.path(hostName,"js","html2canvas.js"),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
		   cat('<script type="text/javascript" src="',file.path(hostName,"js","jspdf","libs","sprintf.js"),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
		   cat('<script type="text/javascript" src="',file.path(hostName,"js","jspdf.js"),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
		   cat('<script type="text/javascript" src="',file.path(hostName,"js","jspdf","libs","base64.js"),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
		   cat('<script type="text/javascript" src="',file.path(hostName,"js","search.js"),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
		   
		   if(!is.null(geneSetNet)){
		  	 cat('<script type="text/javascript" src="',file.path(hostName,"js","cytoscapeweb","js","min","json2.min.js"),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
		   	 cat('<script type="text/javascript" src="',file.path(hostName,"js","cytoscapeweb","js","min","AC_OETags.min.js"),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
			   cat('<script type="text/javascript" src="',file.path(hostName,"js","cytoscapeweb","js","min","cytoscapeweb.min.js"),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
			   cat('<script type="text/javascript" src="',file.path(hostName,"js","cytoscapeWeb.js"),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
			 }
			 
			 if(!is.null(geneSetDAG)){
						 			
			 	 cat('<script type="text/javascript" src="',file.path(hostName,"js","d3","d3.v3.min.js"),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
			 	 cat('<script type="text/javascript" src="',file.path(hostName,"js","bootstrap.min.js"),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")	 
			 	 cat('<script type="text/javascript" src="',file.path(hostName,"js","jquery.tipsy.js"),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
			 	 cat('<script type="text/javascript" src="',file.path(hostName,"js","jquery.contextMenu.js"),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
			 	 cat('<script type="text/javascript" src="',file.path(hostName,"js","jquery.floatThead.min.js"),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
				 cat('<script type="text/javascript" src="',file.path(hostName,"js","dagre","dagre.js"),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
			 	 cat('<script type="text/javascript" src="',file.path(hostName,"js","Minimap.js"),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
				 cat('<script type="text/javascript" src="',file.path(hostName,"js","MinimapZoom.js"),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
				 cat('<script type="text/javascript" src="',file.path(hostName,"js","DirectedAcyclicGraph.js"),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
				 cat('<script type="text/javascript" src="',file.path(hostName,"js","List.js"),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
				 cat('<script type="text/javascript" src="',file.path(hostName,"js","Selectable_svv.js"),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
				 cat('<script type="text/javascript" src="',file.path(hostName,"js","Graph.js"),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
				 cat('<script type="text/javascript" src="',file.path(hostName,"js","History.js"),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
				 cat('<script type="text/javascript" src="',file.path(hostName,"js","Tooltip.js"),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
			   cat('<script type="text/javascript" src="',file.path(hostName,"js","FileSaver.js"),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
			   cat('<script type="text/javascript" src="',file.path(hostName,"js","download.js"),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
			   cat('<script type="text/javascript" src="',file.path(hostName,"js","ContextMenu.js"),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
			   cat('<script type="text/javascript" src="',file.path(hostName,"js","searchView.js"),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
			   cat('<script type="text/javascript" src="',file.path(hostName,"js","canvg","rgbcolor.js"),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
			   cat('<script type="text/javascript" src="',file.path(hostName,"js","canvg","StackBlur.js"),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
			   cat('<script type="text/javascript" src="',file.path(hostName,"js","canvg","canvg.js"),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
			   cat('<script type="text/javascript" src="',file.path(hostName,"js", "xtrace_utils.js"),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
			   cat('<script type="text/javascript" src="',file.path(hostName,"js", "xtrace_graph.js"),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
			   cat('<script type="text/javascript" src="',file.path(hostName,"js", "plotDAG.js"),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
			 }
		}
	 
	  cat('<script type="text/javascript" src="',file.path(hostName,"js","tooltipster.bundle.min.js"),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
	  cat("<script>$(document).ready(function(){$('.tooltip-container').tooltipster({functionInit: function (instance, helper) {var content = $(helper.origin).find('.tooltip_content').detach();instance.content(content);}});});</script>\n",file=outputHtmlFile,append=TRUE,sep="")

    cat('<link href="',file.path(hostName,"css","jquery-ui.css"),'" rel="stylesheet" type="text/css">\n',file=outputHtmlFile,append=TRUE,sep="")
   
    if(!is.null(geneSetDAG)){
	    cat('<link href="',file.path(hostName,"css","xtrace.css"),'" rel="stylesheet" type="text/css">\n',file=outputHtmlFile,append=TRUE,sep="")
	    cat('<link href="',file.path(hostName,"css","tipsy.css"),'" rel="stylesheet" type="text/css">\n',file=outputHtmlFile,append=TRUE,sep="")
	    cat('<link href="',file.path(hostName,"css","jquery.contextMenu.css"),'" rel="stylesheet" type="text/css">\n',file=outputHtmlFile,append=TRUE,sep="")
	    cat('<link href="http://maxcdn.bootstrapcdn.com/bootstrap/3.3.5/css/bootstrap-theme.min.css" rel="stylesheet" type="text/css">\n',file=outputHtmlFile,append=TRUE,sep="")
	    cat('<link href="http://maxcdn.bootstrapcdn.com/bootstrap/3.3.5/css/bootstrap.min.css" rel="stylesheet" type="text/css">\n',file=outputHtmlFile,append=TRUE,sep="")
	 	  cat('<link href="',file.path(hostName,"css","graph.css"),'" rel="stylesheet" type="text/css">\n',file=outputHtmlFile,append=TRUE,sep="")
    }
    cat('<link href="',file.path(hostName,"css","report.css"),'" rel="stylesheet" type="text/css">\n',file=outputHtmlFile,append=TRUE,sep="")
	  cat('<link href="',file.path(hostName,"css","zebra.css"),'" rel="stylesheet" type="text/css">\n',file=outputHtmlFile,append=TRUE,sep="")
	  cat('<link href="',file.path(hostName,"css","search.css"),'" rel="stylesheet" type="text/css">\n',file=outputHtmlFile,append=TRUE,sep="")
	  cat('<link href="',file.path(hostName,"css","tooltipster.bundle.min.css"),'" rel="stylesheet" type="text/css">\n',file=outputHtmlFile,append=TRUE,sep="")

	  cat("</head>\n",file=outputHtmlFile,append=TRUE)
	 
	  cat("<body>\n",file=outputHtmlFile,append=TRUE)
	  cat('<p align="center">\n',file=outputHtmlFile,append=TRUE)
	  cat('<iframe src="',file.path(hostName,"html","head.html"),'" style="border:0px;height:180px;width:100%"></iframe>\n',file=outputHtmlFile,append=TRUE,sep="")
}
