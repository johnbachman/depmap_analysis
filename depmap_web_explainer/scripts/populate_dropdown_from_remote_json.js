$(function(){
    // Load everything

    // Populate first dropdown from json array at:
    // https://s3.amazonaws.com/depmap-public/explainable_ids_1534216288.json
    var first_select_list = "https://s3.amazonaws.com/depmap-public/explainable_ids_1534216288.json";
    var select_first_gene, $select_first_gene
    var select_second_gene, $select_second_gene
    var indra_server_addr = "https://l3zhe2uu9c.execute-api.us-east-1.amazonaws.com/dev/statements/from_hashes";
    var indra_english_asmb = "http://api.indra.bio:8000/assemblers/english";

    function isNumeric(n) {
        return !isNaN(parseFloat(n)) && isFinite(n);
    }

    function isInt(value) {
      return !isNaN(value) && 
             parseInt(Number(value)) == value && 
             !isNaN(parseInt(value, 10));
    }

    function getEnglishByJson(json_stmt_array) {
        eng_stmt = $.ajax({
            url: indra_english_asmb,
            type: "POST",
            dataType: "json",
            contentType: "application/json",
            data: JSON.stringify(json_stmt_array),
        });
        return eng_stmt
    };

    function getStatementByHash(indra_query) {
        stmts_db = $.ajax({
            url: indra_server_addr,
            type: "POST",
            dataType: "json",
            contentType: "application/json",
            data: JSON.stringify(indra_query),
            });
        return stmts_db;
    };

    function grabJSON (url, callback) {
        return $.ajax({url: url, dataType: "json"});
    };

    $select_second_gene = $("#select_second_gene").selectize({
        valueField: "second_item",
        labelField: "second_item",
        searchField: ["second_item"],

        // A single field or an array of fields to sort by.
        sortField: {
            field: "second_item",
            direction: "asc" 
        },

        onChange: function(value) {
                    
            // Refer to the div (or other object) where the output text should be
            let output_text = $("#my_outputB")[0];

            // Add an empty innerHTML object (otherwise it keeps appending to current HTML object)
            output_text.innerHTML = null;
            
            // Build an element (here: "thingy") in the innerHTML that includes the selected value.
            thingy = document.createElement("span");

            // Get the text from the selected item dropdown 
            // thingy.textContent = "Subject: " + subj_input.options[subj_input.selectedIndex].text
            thingy.textContent = "Gene B: " + value

            // Append the element to the div object
            output_text.appendChild(thingy)

            // SET ADDRESSES TO AWS S3 DATA
            // Query of evidence for A->B
            // s3 bucket prefix string
            // Examples:
            // https://s3.amazonaws.com/depmap-public/indra_db_20180730_hash_json/CSDE1_is_subj.json <-- lookup objects given subject
            // https://s3.amazonaws.com/depmap-public/indra_db_20180730_hash_json/CSDE1_is_obj.json <-- lookup subjects given object
            // https://s3.amazonaws.com/depmap-public/correlation_pairs_above_03/correlates_with_A1BG.json
            // https://s3.amazonaws.com/depmap-public/Q3_depmap_20180730_db_explained_improved/A1BG_is_subj.json

            let s3_prefix = "https://s3.amazonaws.com/depmap-public/";
            let s3_subj_expl = "Q3_depmap_20180730_db_explained_improved/";
            let s3_indra_db = "indra_db_20180730_hash_json/"; // INDRA DB LOOKUP
            let s3_correlations = "correlation_pairs_above_03/correlates_with_";

            // Get the current selections
            var geneA = document.getElementById("select_first_gene").value;
            var geneB = document.getElementById("select_second_gene").value;

            // Set adresses
            var geneA_is_subj_expl_address = s3_prefix + s3_subj_expl + geneA + "_is_subj.json";
            var geneB_is_subj_expl_address = s3_prefix + s3_subj_expl + geneB + "_is_subj.json";
            var geneA_is_subj_address = s3_prefix + s3_indra_db + geneA + "_is_subj.json";
            var geneB_is_subj_address = s3_prefix + s3_indra_db + geneB + "_is_subj.json";
            var geneA_is_obj_address = s3_prefix + s3_indra_db + geneA + "_is_obj.json";
            var geneB_is_obj_address = s3_prefix + s3_indra_db + geneB + "_is_obj.json";
            var correlates_with_A = s3_prefix + s3_correlations + geneA + ".json";
            var depmap1 = "https://depmap.org/portal/interactive/?xDataset=Avana&xFeature="
            var depmap2 = "&yDataset=Avana&yFeature="
            var depmap3 = "&colorDataset=lineage&colorFeature=all&filterDataset=context&filterFeature=&regressionLine=false&statisticsTable=false&associationTable=true&plotOnly=false"

            // Query and output A, B, correlation
            var A_B_correlation = $.ajax({
                url: correlates_with_A,
                success: function(res) {
                    correlation_AB = res[geneB]

                    var correlation_output = $("#show_correlation")[0];
                    correlation_output.innerHTML = null;
                    var correlation_output_element = document.createElement("a")
                    var linkText = document.createTextNode("Link to depmap plot")
                    correlation_output_element.appendChild(linkText);
                    correlation_output_element.title = "Link to depmap plot for " + geneA + " vs " + geneB
                    correlation_output_element.href = depmap1 + geneA + depmap2 + geneB + depmap3
                    if (isNumeric(correlation_AB)) {
                        correlation_output_element.textContent = geneA + ", " + geneB + ", " + parseFloat(correlation_AB).toFixed(4).toString() // DECIMAL PLACES IN CORRELATION
                    } else {
                        // When we don't have the correation; If it happens, you probably need to update the correlation jsons
                        console.log('Correlation is not a valid number!')
                        correlation_output_element.textContent = geneA + ", " + geneB + ", (not available)"
                    }
                    correlation_output.appendChild(correlation_output_element)
                },
                error: function() {
                    var correlation_output = $("#show_correlation")[0];
                    correlation_output.innerHTML = null;
                    var correlation_output_element = document.createElement("span");
                    correlation_output_element.textContent = "Failed to load from " + correlates_with_A
                    correlation_output.appendChild(correlation_output_element)
                }
            })

            // To be used so we can query common up/downstream on B-X-A when A->B gives back a result but not B->A;
            // Should also use it for avoiding double output.
            var AcB_d_output = false
            var AB_im_output = false
            // var BA_im_output = false

            // Query and output all subj:A -> obj:B
            var geneA_is_subj_promise = $.ajax({
                url: geneA_is_subj_expl_address,
                success: function(res) {
                    let obj = geneB

                    // Should return a dict of the format below
                    connection_type_list = res[obj]
                    
                    // json Return format:
                    // {"CHEK1":
                    //          {'directed': [["Phosphorylation", 17052011326019041], ["Phosphorylation", -32662422560218481],
                    //                        ["Dephosphorylation", -750973640471511], ["Activation", 30186062639888508],
                    //                        ["Inhibition", 20888016729018787], ["Activation", -8364720323695997]],
                    //           'undirected': ["Complex", 35575713738557636], 
                    //           'x_is_intermediary': [X],
                    //           'x_is_downstream': [X],
                    //           'x_is_upstream': [X]}}
                    //           
                    // access like:
                    // responseJSON["HGNC_id"][n][0/1]
                    // where n = number of interaction types (7 in above example) and [0/1] will give
                    // "type"/statement hash.

                    // OUTPUT EXPLANATIONS

                    // if connection undirected
                    if (connection_type_list.undirected.length > 0) {
                        var debug_string = 'output_AcB' // Kept for now in anticipation of future debugging needs
                        // console.log(debug_string)

                        AcB_d_output = true

                        // Reference and initialize the output pointer
                        var output_AcB = $("#expl_A_complex_B")[0];
                        output_AcB.innerHTML = null;

                        // Get reference to the text badge so we can output evidence count
                        var AcB_ev_count = document.getElementById("collapseAcB_ev_count");

                        // args :      output_pointer, ev_counter_pointer, type_hash_array, debug_string
                        output_directs(output_AcB, AcB_ev_count, connection_type_list.undirected, debug_string)
                    }

                    // if connection directed
                    if (connection_type_list.directed.length > 0) {
                        var debug_string = 'output_AB'
                        // console.log(debug_string)

                        // Reference and initialize the output pointer
                        var output_AB = $("#expl_A_to_B")[0];
                        output_AB.innerHTML = null;

                        // Get reference to the text badge so we can output evidence count
                        var AB_ev_count = document.getElementById("collapseAB_ev_count");

                        output_directs(output_AB, AB_ev_count, connection_type_list.directed, debug_string)
                    }

                    // 'x_is_intermediary'; This is for A->X->B; B->X->A is/should be below in 'geneB_is_subj_promise'
                    if (connection_type_list.x_is_intermediary.length > 0) {
                        var debug_string = 'output_AXB'
                        // console.log(debug_string)

                        var output_AXB = $("#expl_A_to_X_to_B")[0];
                        output_AXB.innerHTML = null;

                        // Get pointer to evidence counter
                        var AXB_ev_count = document.getElementById("collapseAXB_ev_count");

                        // output_intermediary(output_pointer, ev_counter_pointer, X_array, subj, obj, SUBJ_is_subj_address, OBJ_is_obj_address, debug_string)
                        output_intermediary(output_AXB, AXB_ev_count, connection_type_list.x_is_intermediary, geneA, geneB, geneA_is_subj_address, geneB_is_obj_address, debug_string)
                    }

                    // 'x_is_downstream'
                    if (connection_type_list.x_is_downstream.length > 0) {
                        var debug_string = 'output_ABx'
                        // console.log(debug_string)

                        AB_im_output = true
                        var output_ABx = $('#expl_x_is_downstream')[0];
                        output_ABx.innerHTML = null;

                        // evidence count pointer
                        var ABx_ev_count = document.getElementById("collapse_st_X_count");

                        // args : output_pointer, ev_counter_pointer, X_array, geneA, geneB, geneA_is_subj_address, geneB_is_subj_address, debug_string
                        output_x_is_downstream(output_ABx, ABx_ev_count, connection_type_list.x_is_downstream, geneA, geneB, geneA_is_subj_address, geneB_is_subj_address, debug_string)
                    }

                    // 'x_is_upstream':
                    if (connection_type_list.x_is_upstream.length > 0) {
                        var debug_string = 'output_xAB'
                        // console.log(debug_string)

                        AB_im_output = true
                        var output_xAB = $('#expl_x_is_upstream')[0];
                        output_xAB.innerHTML = null;

                        // evidence count pointer
                        var xAB_ev_count = document.getElementById("collapse_sr_X_count");

                        // args : output_pointer, ev_counter_pointer, X_array, geneA, geneB, geneA_is_obj_address, geneB_is_obj_address, debug_string
                        output_x_is_upstream(output_xAB, xAB_ev_count, connection_type_list.x_is_upstream, geneA, geneB, geneA_is_obj_address, geneB_is_obj_address, debug_string)
                    }
                },
                error: function() {
                    var output_AB = $("#expl_A_to_B")[0];
                    output_AB.innerHTML = null;
                    let AB_output_element_err = document.createElement("div")
                    AB_output_element_err.textContent = "Could not query " + geneA_is_subj_expl_address
                    output_AB.appendChild(AB_output_element_err)
                }

            })

            // Query and output all subj:B -> obj:A
            var geneB_is_subj_promise = $.ajax({
                url: geneB_is_subj_expl_address,
                success: function(res) {
                    let obj = geneA
                    connection_type_list = res[obj]

                    // if connection undirected and not already printed
                    if (!AcB_d_output) {
                        if (connection_type_list.undirected.length > 0) {
                            var debug_string = 'output_AcB'
                            // console.log(debug_string)

                            AcB_d_output = true

                            // Reference and initialize the output pointer
                            var output_AcB = $("#expl_A_complex_B")[0];
                            output_AcB.innerHTML = null;

                            // Get reference to the text badge so we can output evidence count
                            var AcB_ev_count = document.getElementById("collapseAcB_ev_count");

                            output_directs(output_AcB, AcB_ev_count, connection_type_list.undirected, debug_string)
                        }
                    }

                    // if connection directed
                    if (connection_type_list.directed.length > 0) {
                        var debug_string = 'output_BA'
                        // console.log(debug_string)

                        var output_BA = $("#expl_B_to_A")[0];
                        output_BA.innerHTML = null;

                        // Evidence counter
                        collapseAB_ev_count
                        var BA_ev_count = document.getElementById("collapseBA_ev_count");

                        // args : output_pointer, ev_counter_pointer, type_hash_array, debug_string
                        output_directs(output_BA, BA_ev_count, connection_type_list.directed, debug_string)
                    }

                    // 'x_is_intermediary'; B->X->A
                    if (connection_type_list.x_is_intermediary.length > 0) {
                        var debug_string = 'output_BXA'
                        // console.log(debug_string)

                        var output_BXA = $("#expl_B_to_X_to_A")[0];
                        output_BXA.innerHTML = null;

                        // Get pointer to evidence counter
                        var BXA_ev_count = document.getElementById("collapseBXA_ev_count");

                        // output_intermediary(output_pointer, ev_counter_pointer, X_array, subj, obj, SUBJ_is_subj_address, OBJ_is_obj_address, debug_string)
                        output_intermediary(output_BXA, BXA_ev_count, connection_type_list.x_is_intermediary, geneB, geneA, geneB_is_subj_address, geneA_is_obj_address, debug_string)
                    }

                    // Check if any output already is up for xAB or ABx
                    if (!AB_im_output) {

                        // 'x_is_downstream'
                        if (connection_type_list.x_is_downstream.length > 0) {
                            var debug_string = 'output_BAx'
                            // console.log(debug_string)

                            var output_ABx = $('#expl_x_is_downstream')[0];
                            output_ABx.innerHTML = null;

                            // evidence count pointer
                            var ABx_ev_count = document.getElementById("collapse_st_X_count");

                            // args : output_pointer, ev_counter_pointer, X_array, geneA, geneB, geneA_is_subj_address, geneB_is_subj_address, debug_string
                            output_x_is_downstream(output_ABx, ABx_ev_count, connection_type_list.x_is_downstream, geneA, geneB, geneA_is_subj_address, geneB_is_subj_address, debug_string)
                        }

                        // 'x_is_upstream':
                        if (connection_type_list.x_is_upstream.length > 0) {
                            var debug_string = 'output_xBA'
                            // console.log(debug_string)

                            var output_xAB = $('#expl_x_is_upstream')[0];
                            output_xAB.innerHTML = null;

                            // evidence count pointer
                            var xAB_ev_count = document.getElementById("collapse_sr_X_count");

                            // args : output_pointer, ev_counter_pointer, X_array, geneA, geneB, geneA_is_obj_address, geneB_is_obj_address, debug_string
                            output_x_is_upstream(output_xAB, xAB_ev_count, connection_type_list.x_is_upstream, geneA, geneB, geneA_is_obj_address, geneB_is_obj_address, debug_string)
                        }
                    }
                },
                error: function() {
                    var output_BA = $("#expl_B_to_A")[0];
                    output_BA.innerHTML = null;
                    let BA_output_element_err = document.createElement("div")
                    BA_output_element_err.textContent = "Could not query " + geneB_is_subj_expl_address
                    output_BA.appendChild(BA_output_element_err)
                }
            })

        } // Closing bracket for onChange: function().
          // Some things that are not defined outside of here: geneA, geneB, all the query addresses
    }); // Closing bracket for $("#select_second_gene").selectize()

    // Save response in a promise
    var json_promise = grabJSON(first_select_list);

    // Use the .then() of promise to make the call only after we have a response.
    // By some JS magic, res == json_promise.responseJSON instead of res == json_promise
    json_promise.then(function(res) {
        data = res
        var items = data.map(function(x) { return { item: x }; })

        $select_first_gene = $("#select_first_gene").selectize({
            // Selctizre usage documentation
            // https://github.com/selectize/selectize.js/blob/master/docs/usage.md

            // Populate the dropdown options from an array ("items")
            options: items,

            // The name of the property in "items" to populate the dropdown with.
            labelField: "item",

            // The name of the property to use as the value when an item is selected.
            valueField: "item",

            // Allows the user to create new items not in the initial list.
            // create: true,

            // An array of property names to analyze when filtering options.
            searchField: ["item"],

            // A single field or an array of fields to sort by.
            sortField: {
                field: "item",
                direction: "asc" 
            },

            // dropdownParent: "body",

            // Updates the current selection 
            onChange: function(value) {
                // Refer to the div (or other object) where the output text should be
                let output_text = $("#my_outputA")[0];

                // Add an empty innerHTML object (otherwise it keeps appending to current HTML object)
                output_text.innerHTML = null;
                
                // build an element (here: "thingy") in the innerHTML that includes the selected value.
                let thingy = document.createElement("span");

                // Get the text from the selected item dropdown 
                thingy.textContent = "Gene A: " + value

                // Append the element to the div object
                output_text.appendChild(thingy)
                
                // Set second query address example:
                // https://s3.amazonaws.com/depmap-public/prior_filtered_neighbor_lookup/neighbors_to_BRCA1.json
                // https://s3.amazonaws.com/depmap-public/neighbor_lookup/neighbors_to_A1BG.json
                s3_prefix = "https://s3.amazonaws.com/depmap-public/neighbor_lookup/neighbors_to_";
                var second_dd_address = s3_prefix + value + ".json"
                // console.log(second_dd_address)

                // Query for next dropdown
                if (!value.length) return;
                select_second_gene.disable();
                select_second_gene.clearOptions();
                select_second_gene.load(function(callback) {
                    var second_json = $.ajax({
                        url: second_dd_address,
                        success: function(results) {
                            var second_items = results.map(function(x) { return {second_item: x }; })
                            select_second_gene.enable();
                            callback(second_items);
                        },
                        error: function() {
                            let output_text = $("#my_outputB")[0];
                            output_text.innerHTML = null;
                            let output_text_err = document.createElement("div")
                            output_text_err.textContent = "Could not load from " + second_dd_address
                            output_text.appendChild(output_text)
                        }
                    })
                });
            } // This is closing bracket for "onChange: function(value)"
        }); // This is closing bracket for "$("#select_first_gene").selectize"

        var select_first_gene = $select_first_gene[0].selectize;
        var select_second_gene = $select_second_gene[0].selectize;

    });

    // Function for quering and outputting plain english description and statement evidence with PMIDs
    function output_directs(output_pointer, ev_counter_pointer, type_hash_array, debug_string){

        // Create array to store each statement hash
        let hash_list = [];

        // debugger; // type_hash_array

        // type_hash_array contains [["type1", hash1], ["type2", hash2], ...]
        for (let i = 0; i < type_hash_array.length; i++) {
            hash_list.push(type_hash_array[i][1]);
        }

        let hash_query = {"hashes": hash_list}
        // let stmt_promise = getStatementByHash(hash_query)
        var stmt_promise = getStatementByHash(hash_query)
        
        stmt_promise.then(function(stmt_response){
            // Get statements array
            var stmts = stmt_response.statements

            // Check error for 'Object.keys(stmts).map(function(x) {return Number(x)})': 
            // jQuery.Deferred exception: Cannot convert undefined or null to object TypeError: Cannot convert undefined or null to object at Function.keys
            // debugger;

            // Generating new hash list here in case the respones are not exactly corresponding to the hash query
            var hash_list = Object.keys(stmts)

            // We could send an array of statement jsons, but then we 
            // would have to keep track of which uuid is with which statement
            // because I don't know if they're being returned in the same order
            // as they were sent in. Instead, let's loop over statements and
            // IDs for now

            // Array to store query responses and corresponding uuids
            var eng_res_array = [];
            var stmt_uuid_array = [];

            // Loop hashes for stmt jsons and store uuid and plain english query response
            for (let hash of hash_list) {
                stmt_json = stmts[hash]
                stmt_uuid_array.push(stmt_json.id)
                json_stmt_array = {"statements": [stmt_json]}
                eng_res_array.push(getEnglishByJson(json_stmt_array))
            }

            // Array Promises; For docs, see:
            // https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Promise/all
            Promise.all(eng_res_array).then(function(eng_array) {
                // Output statement count
                let output_element_stmt_count = document.createElement("h4")
                output_element_stmt_count.textContent = "Found " + eng_array.length + " statements."
                output_element_stmt_count.style = "background-color:#F2F2F2;"
                output_pointer.appendChild(output_element_stmt_count)

                // Update the count in the badge if it is a 'real' pointer (the A-X-B won't have evidence counter for now)
                ev_counter_pointer.textContent = eng_array.length  // EVIDENCE SOURCE COUNT
                
                // Loop response array for plain english
                for (let k = 0; k < eng_array.length; k++) {
                    // Get uuid, english output
                    sid = stmt_uuid_array[k]
                    eng_res = eng_array[k]
                    eng_plain = eng_res.sentences[sid]

                    // Count evidence
                    stmt_json = stmts[hash_list[k]]
                    ev_len = stmt_json.evidence.length

                    // Output for Plain English
                    let output_element_pe = document.createElement("h4")
                    // output_element_pe.style = "display:inline-block; margin-right:10px;"; // For placement of text and buttons
                    output_element_pe.textContent = (k+1) + ": " + eng_plain + ' (' + ev_len + ' sources)'
                    output_pointer.appendChild(output_element_pe)

                    // EVIDENCE BUTTON
                    let ev_button_div = document.createElement("div");
                    // ev_button_div.style = "display:inline-block; margin-right:10px;";
                    let ev_button_output_text = document.createElement("div");
                    ev_button_output_text.id = sid;
                    // ev_button_output_text.style = "display:inline-block; margin-right: 10px;";
                    ev_button_output_text.textContent = "Click on the button far to the right :) "
                    let ev_button = document.createElement("button");
                    ev_button.classList.add("btn", "btn-default", "btn-evidence", "pull-right");
                    ev_button.textContent = "clickme"; // BUTTON TEXT
                    ev_button.dataset.index = k // Store index
                    ev_button.dataset.id = sid; // Button identifier
                    ev_button_div.appendChild(ev_button)
                    ev_button_div.appendChild(ev_button_output_text)
                    output_pointer.appendChild(ev_button_div)

                }

                $(".btn-evidence").on("click", function(b){
                    // Loop through the evidence for the statement the button is linked to
                    n = b.currentTarget.dataset.index
                    btn_id = b.currentTarget.dataset.id
                    hash = hash_list[n]
                    var stmt_json = stmts[hash]
                    var ev_output_div = $("#"+btn_id)[0];
                    ev_output_div.innerHTML = null;

                    for (let k = 0; k < stmt_json.evidence.length; k++) {

                        _pmid = stmt_json.evidence[k].pmid

                        // Output for pmid link
                        let output_element_pmid = document.createElement("a")
                        output_element_pmid.href = "https://www.ncbi.nlm.nih.gov/pubmed/" + _pmid
                        output_element_pmid.textContent = "PMID: " + _pmid
                        ev_output_div.appendChild(output_element_pmid)

                        // Ouput for all evidence of of current stmt
                        let output_element_ev = document.createElement("div")
                        output_element_ev.innerHTML = null;
                        output_element_ev.textContent = "\"" + stmt_json.evidence[k].text + "\""
                        ev_output_div.appendChild(output_element_ev)
                    }
                });
            });
        })
    } // Closes the output_directs function bracket

    // Use this function for s->X->o (same for both ways, so make two calls) the query needs to be over two jsons: SUBJ_is_subj and OBJ_is_obj
    function output_intermediary(output_pointer, x_counter_pointer, x_array, subj, obj, SUBJ_is_subj_address, OBJ_is_obj_address, debug_string){
        var rand_id = Number(Math.random()*10**17).toString(); // Just create a random id
        dropdown_div = document.createElement("div");
        dropdown_div.id = rand_id;
        dropdown_div.style = "width: 520px; top: 36px; left: 0px; visibility: visible;"
        output_pointer.appendChild(dropdown_div)
        var items = x_array.map(function(x) { return { item: x }; })

        // Update the count of X in the badge
        cur_count = x_counter_pointer.textContent
        if (isInt(cur_count)) {
            x_counter_pointer.textContent = x_array.length
        } else {
            console.log('Could not interpret badge count as integer!')
            x_counter_pointer.textContent = cur_count
        }

        // Create dropdown with all X
        $select_intermediate = $("#"+rand_id).selectize({
            options: items,
            valueField: "item",
            labelField: "item",
            searchField: ["item"],

            // A single field or an array of fields to sort by.
            sortField: {
                field: "item",
                direction: "asc" 
            },

            // On select: Query A-X and B-X and output the english statements and their evidence
            onChange: function(x_value) {
                if (!x_value.length) return;
                two_promises = [];
                two_promises.push(grabJSON(SUBJ_is_subj_address)) // <--- Query for 'SUBJ_is_subj'
                two_promises.push(grabJSON(OBJ_is_obj_address)) // <--- Query for 'OBJ_is_obj'

                // Wait for both promises to resolve
                Promise.all(two_promises).then(function(two_jsons_ar){
                    // Get the the hash arrays
                    S_is_subj_lookup = two_jsons_ar[0];
                    O_is_obj_lookup = two_jsons_ar[1];

                    // args : output_pointer, ev_counter_pointer, type_hash_array
                    output_directs(output_pointer, x_counter_pointer, S_is_subj_lookup[x_value], debug_string)
                    output_directs(output_pointer, x_counter_pointer, O_is_obj_lookup[x_value], debug_string)

                });
            }
        })

    } // Closes the output_intermediary function bracket

    // Here we need 'A_is_obj' and 'B_is_obj'
    function output_x_is_upstream(output_pointer, x_counter_pointer, x_array, geneA, geneB, geneA_is_obj_address, geneB_is_obj_address, debug_string){
        var rand_id = Number(Math.random()*10**17).toString(); // Just create a random id
        dropdown_div = document.createElement("div");
        dropdown_div.id = rand_id;
        dropdown_div.style = "width: 520px; top: 36px; left: 0px; visibility: visible;"
        output_pointer.appendChild(dropdown_div)
        var items = x_array.map(function(x) { return { item: x }; })

        // Update the count of X in the badge
        cur_count = x_counter_pointer.textContent
        if (isInt(cur_count)) {
            x_counter_pointer.textContent = x_array.length
        } else {
            console.log('Could not interpret badge count as integer!')
            x_counter_pointer.textContent = cur_count
        }

        $select_upstream_x = $("#"+rand_id).selectize({
            options: items,
            valueField: "item",
            labelField: "item",
            searchField: ["item"],

            // A single field or an array of fields to sort by.
            sortField: {
                field: "item",
                direction: "asc" 
            },

            // On select: Query A-X and B-X and output the english statements and their evidence
            onChange: function(x_value) {
                if (!x_value.length) return;
                // Query A-X and B-X
                two_promises = [];
                two_promises.push(grabJSON(geneA_is_obj_address))
                two_promises.push(grabJSON(geneB_is_obj_address))

                // Wait for both promises to resolve
                Promise.all(two_promises).then(function(two_jsons_ar){
                    A_is_obj_lookup = two_jsons_ar[0];
                    B_is_obj_lookup = two_jsons_ar[1];

                    // args : output_pointer, ev_counter_pointer, type_hash_array
                    output_directs(output_pointer, x_counter_pointer, A_is_obj_lookup[x_value], debug_string);
                    output_directs(output_pointer, x_counter_pointer, B_is_obj_lookup[x_value], debug_string);
                });
            }
        })
    }

    // Here we need 'A_is_subj' and 'B_is_subj'
    function output_x_is_downstream(output_pointer, x_counter_pointer, x_array, geneA, geneB, geneA_is_subj_address, geneB_is_subj_address, debug_string){
        var rand_id = Number(Math.random()*10**17).toString(); // Create a random id
        dropdown_div = document.createElement("div");
        dropdown_div.id = rand_id;
        dropdown_div.style = "width: 520px; top: 36px; left: 0px; visibility: visible;"
        output_pointer.appendChild(dropdown_div)
        var items = x_array.map(function(x) { return { item: x }; })

        // Update the count of X in the badge
        cur_count = x_counter_pointer.textContent
        if (isInt(cur_count)) {
            x_counter_pointer.textContent = x_array.length
        } else {
            console.log('Could not interpret badge count as integer!')
            x_counter_pointer.textContent = cur_count
        }

        $select_downstream_x = $("#"+rand_id).selectize({
            options: items,
            valueField: "item",
            labelField: "item",
            searchField: ["item"],

            // A single field or an array of fields to sort by.
            sortField: {
                field: "item",
                direction: "asc" 
            },

            // On select: Query A-X and B-X and output the english statements and their evidence
            onChange: function(x_value) {
                if (!x_value.length) return;
                // Query A-X and B-X
                two_promises = [];
                two_promises.push(grabJSON(geneA_is_subj_address));
                two_promises.push(grabJSON(geneB_is_subj_address));

                // Wait for both promises to resolve
                Promise.all(two_promises).then(function(two_jsons_ar){
                    A_is_subj_lookup = two_jsons_ar[0];
                    B_is_subj_lookup = two_jsons_ar[1];

                    // args : output_pointer, ev_counter_pointer, type_hash_array, debug_string
                    output_directs(output_pointer, x_counter_pointer, A_is_subj_lookup[x_value], debug_string);
                    output_directs(output_pointer, x_counter_pointer, B_is_subj_lookup[x_value], debug_string);
                });
            }
        })
    }
}); // This closes the load everything bracket
