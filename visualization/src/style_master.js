/* style_master.js
 *
 * Program the visual styles for Omics Integrator visualization
 * using the Cytoscape.js library.
 *
 */

var style = [ {
  "format_version" : "1.0",
  "target_cytoscapejs_version" : "~2.1",
  "title" : "default",
  "style" : [ {
    "selector" : "edge",
    "css" : {
      "target-arrow-shape" : "none",
      "opacity" : 1.0,
      "content" : "",
      "font-size" : 10,
      "text-opacity" : 1.0,
      "target-arrow-color" : "rgb(64,64,64)",
      "line-color" : "rgb(64,64,64)",
      "font-family" : "Dialog.plain",
      "font-weight" : "normal",
      "source-arrow-shape" : "none",
      "color" : "rgb(0,0,0)",
      "source-arrow-color" : "rgb(64,64,64)",
      "width" : 2.0,
      "line-style" : "solid",
      "curve-style": "bezier"
    }
  }, {
    "selector" : "edge[interaction = 'pp']",
    "css" : {
      "target-arrow-shape" : "none"
    }
  }, {
    "selector" : "edge[interaction = 'pd']",
    "css" : {
      "target-arrow-shape" : "triangle",
      "target-arrow-fill" : "filled"
    }
  }, {
    "selector" : "edge[interaction = 'tp']",
    "css" : {
      "target-arrow-shape" : "triangle",
      "target-arrow-fill" : "filled"
    }
  }, {
    "selector" : "edge[interaction = 'pp']",
    "css" : {
      "line-color" : "rgb(153,153,153)",
      "target-arrow-color" : "rgb(153,153,153)",
      "source-arrow-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'pd']",
    "css" : {
      "line-color" : "rgb(51,153,0)",
      "target-arrow-color" : "rgb(51,153,0)",
      "source-arrow-color" : "rgb(51,153,0)"
    }
  }, {
    "selector" : "edge[interaction = 'tp']",
    "css" : {
      "line-color" : "rgb(152,233,7)",
      "target-arrow-color" : "rgb(152,233,7)",
      "source-arrow-color" : "rgb(152,233,7)"
    }
  },{
    "selector" : "edge:selected",
    "css" : {
      "line-color" : "rgb(255,0,0)"
    }
  }, {
    "selector" : "node",
    "css" : {
      "color" : "rgb(0,0,0)",
      "background-opacity" : 1.0,
      "border-color" : "gray",
      "text-opacity" : 1.0,
      "background-color" : "rgb(221,255,255)",
      "width" : 75.0,
      "font-family" : "SansSerif.plain",
      "font-weight" : "normal",
      "text-valign" : "center",
      "text-halign" : "center",
      "border-opacity" : 1.0,
      "border-width" : "2",
      "shape" : "hexagon",
      "height" : 35.0,
      "font-size" : 12,
      "content" : "data(name)"
    }
  }, {
    "selector" : "node[TerminalType = 'TF']",
    "css" : {
      "shape" : "triangle"
    }
  }, {
    "selector" : "node[TerminalType = 'TF']",
    "css" : {
      "width" : 85.0
    }
  }, {
    "selector" : "node[TerminalType = 'mRNA']",
    "css" : {
      "shape" : "diamond"
    }
  }, {
    "selector" : "node[TerminalType = 'mRNA']",
    "css" : {
      "width" : 85.0,
      "height" : 40.0
    }
  }, {
    "selector" : "node[TerminalType = 'Metabolite']",
    "css" : {
      "shape" : "ellipse"
    }
  }, {
    "selector" : "node[TerminalType = 'Metabolite']",
    "css" : {
      "width" : 82.0
    }
  }, {
    "selector" : "node[TerminalType = 'Proteomic']",
    "css" : {
      "shape" : "rectangle"
    }
  }, {
    "selector" : "node[TerminalType = 'mRNA']",
    "css" : {
      "background-color" : "green"
    }
  }, {
    "selector" : "node[ProteinChange < 0]",
    "css" : {
      "background-color" : "blue"
    }
  }, {
    "selector" : "node[ProteinChange > 0]",
    "css" : {
      "background-color" : "red"
    }
  }, {
    "selector" : "node[ProteinChange = 0]",
    "css" : {
      "background-color" : "rgb(221,255,255)"
    }
  }, {
    "selector" : "node[TerminalType = 'Metabolite']",
    "css" : {
      "border-color" : "gray"
    }
  }, {
    "selector" : "node[TerminalType = 'Metabolite']",
    "css" : {
      "border-width" : "2"
    }
  }, {
    "selector" : "node[GeneChange > 0]",
    "css" : {
      "border-color" : "rgb(255,153,0)",
      "border-width" : "7"
    }
  }, {
    "selector" : "node[GeneChange < 0]",
    "css" : {
      "border-color" : "purple",
      "border-width" : "7"
    }
  }, {
    "selector" : "node[GeneChange = 0]",
    "css" : {
      "border-color" : "gray",
      "border-width" : "2"
    }
  }, {
    "selector" : "node[TerminalType = 'mRNA']",
    "css" : {
      "background-blacken" : -0.4
    }
