import json
import sys

def process_graph_json(json_file, graph_filename):
        """
        Formats graph JSON file into javascript code, and extracts min and max values for style

        Args:
                json_file: network JSON filename
                graph_filename: a new file to write the javascript to

        Returns:
                a tuple containing min and max values
                
        """

        #Add the beginning to js file
        graph_file = open(graph_filename, "w")
        graph_file.write("var graph = \n")
        
        #Load JSON file
        json_data = open(json_file)
        network_data = json.load(json_data)
        json_data.close()
	
        node_list = network_data["elements"]["nodes"]
        edge_list = network_data["elements"]["edges"]
	    
        #Iterate through nodes and to keep track and min/max
        #and change clusters to "parent" attribute
        min_pc = 0
        max_pc = 0
        clusters = []
        for i, node in enumerate(node_list):
            cluster = [value for key, value in node["data"].items() if 'cluster' in key.lower()]
            if len(cluster) == 1:
                network_data["elements"]["nodes"][i]["data"]["parent"] = cluster[0]
                if cluster[0] not in clusters:
                    clusters.append(cluster[0])
            if "ProteinChange" in node["data"]:
                current_value = node["data"]["ProteinChange"]
            else:
                continue
            if current_value > max_pc:
                max_pc = current_value
            if current_value < min_pc:
                min_pc = current_value
        
        #Iterate through edges and to keep track and min/max
        min_edge = 0
        max_edge = 0
        for edge in edge_list:
            current_value = edge["data"]["cost"]
            if current_value > max_edge:
                max_edge = current_value
            if current_value < min_edge:
                min_edge = current_value


        #Add parent nodes
        if len(clusters)>0:
            for clus_name in clusters:
                network_data["elements"]["nodes"].append({"data":{"id": clus_name}})

        #Write JSON to js file
        json.dump(network_data, graph_file, indent=4)

        #Add end of js file
        graph_file.write(";")
        graph_file.close()

	
        return [min_pc, max_pc, min_edge, max_edge]

def create_style_code_cytoscape(min_max, style_filename):
        """
        Creates the style JavaScript code to be used for the visualization
        using the cytoscape.js library.

        Args:
                min_max: a tuple containing min and max values
                style_filename: style filename

        Returns:
                None
        """
        
        min_pc, max_pc, min_edge, max_edge = [str(item) for item in min_max]
       
        style_file = open(style_filename, "a")

        style_file.write("\n")

        # ProteinChange
        if max_pc > 0:
            style_file.write("  }, {\n"
                             "    \"selector\" : \"node[ProteinChange > " + max_pc + "]\",\n"
                             "    \"css\" : {\n"
                             "      \"background-blacken\" : 0.5\n"
                             "    }\n"
                             "  }, {\n"
                             "    \"selector\" : \"node[ProteinChange = " + max_pc + "]\",\n"
                             "    \"css\" : {\n"
                             "      \"background-blacken\" : 0.5\n"
                             "    }\n"
                             "  }, {\n"
                             "    \"selector\" : \"node[ProteinChange > 0][ProteinChange < " + max_pc + "]\",\n"
                             "    \"css\" : {\n"
                             "      \"background-blacken\" : \"mapData(ProteinChange, 0, " + max_pc + ", -0.9, 0.5)\"\n"
                             "    }\n")
        if min_pc < 0:
            style_file.write("  }, {\n"
                             "    \"selector\" : \"node[ProteinChange > " + min_pc + "][ProteinChange < 0]\",\n"
                             "    \"css\" : {\n"
                             "      \"background-blacken\" : \"mapData(ProteinChange, " + min_pc + ", 0, 0.5, -0.9)\"\n"
                             "    }\n"
                             "  }, {\n"
                             "    \"selector\" : \"node[ProteinChange = " + min_pc + "]\",\n"
                             "    \"css\" : {\n"
                             "      \"background-blacken\" : 0.5\n"
                             "    }\n"
                             "  }, {\n"
                             "    \"selector\" : \"node[ProteinChange < " + min_pc + "]\",\n"
                             "    \"css\" : {\n"
                             "      \"background-blacken\" : 0.5\n"
                             "    }\n")
        # Edge Width
        style_file.write("  }, {\n"
                             "    \"selector\" : \"edge\",\n"
                             "  \"css\" : {\n"
                             "    \"width\" : \"mapData(cost," + max_edge + "," + min_edge + ",1,10)\"\n"
                             "  }\n")

        # Ending part
        style_file.write("  }, {\n"
                         "    \"selector\" : \"node:selected\",\n"
                         "    \"css\" : {\n"
                         "      \"background-color\" : \"rgb(255,255,0)\"\n"
                         "    }\n"
                         "  } ]\n"
                         "} ];")

        
        style_file.close()        
                      
        return None




if __name__ == '__main__':

        outdir = sys.argv[1]
        user_dir = sys.argv[2]

        min_max_values = process_graph_json("%s/graph_json.json"%user_dir, "%s/graph.js"%outdir)

        create_style_code_cytoscape(min_max_values, "%s/style.js"%outdir)

