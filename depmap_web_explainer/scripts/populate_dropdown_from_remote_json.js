$(function(){
    // Load everything

    // Populate first dropdown from json array at:
    // https://s3.amazonaws.com/depmap-public/explainable_ids.json
    var first_select_list = "https://s3.amazonaws.com/depmap-public/explainable_ids.json";
    var select_first_gene, $select_first_gene
    var select_second_gene, $select_second_gene
    var indra_server_addr = "https://l3zhe2uu9c.execute-api.us-east-1.amazonaws.com/dev/statements/from_hashes";
    var indra_english_asmb = "http://api.indra.bio:8000/assemblers/english";

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

            // Query of evidence for A->B
            // s3 bucket prefix string
            // Examples:
            // https://s3.amazonaws.com/depmap-public/indra_db_20180730_hash_json/41137subj_CSDE1.json
            // https://s3.amazonaws.com/depmap-public/correlation_pairs_above_03/correlates_with_A1BG.json
            // 
            let s3_prefix = "https://s3.amazonaws.com/depmap-public/";
            let s3_subj_interactions = "indra_db_20180730_hash_json/41137subj_";
            let s3_obj_interactions = "indra_db_20180730_hash_json/41137obj_"; // Can be used to show subjects of given object
            let s3_correlations = "correlation_pairs_above_03/correlates_with_";

            // Get the current selections
            var geneA = document.getElementById("select_first_gene").value;
            var geneB = document.getElementById("select_second_gene").value;

            // Set adresses
            var geneA_is_subj_address = s3_prefix + s3_subj_interactions + geneA + ".json";
            var geneB_is_subj_address = s3_prefix + s3_subj_interactions + geneB + ".json";
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
                    correlation_output_element.textContent = geneA + ", " + geneB + ", " + correlation_AB
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

            // Query and output all subj:A -> obj:B
            var geneA_is_subj_promise = $.ajax({
                url: geneA_is_subj_address,
                success: function(res) {
                    let obj = geneB

                    // Should return a nested list of the format above
                    type_hash_list = res[obj]
                    // console.log(type_hash_list)
                    
                    // json Return format:
                    // {"CHEK1":
                    // [["Complex", 35575713738557636], ["Phosphorylation", 17052011326019041], 
                    //  ["Phosphorylation", -32662422560218481], ["Dephosphorylation", -750973640471511],
                    //  ["Activation", 30186062639888508], ["Inhibition", 20888016729018787], 
                    //  ["Activation", -8364720323695997]]}
                    // access like:
                    // responseJSON["HGNC_id"][n][0/1]
                    // where n = number of interaction types (7 in aboverexample) and [0/1] will give
                    // "type"/statement hash.
                    
                    var output_AB = $("#expl_A_to_B")[0];
                    output_AB.innerHTML = null;

                    // Create array to store each statement hash
                    var hash_list = [];
                    var hash_list_len = 0;

                    for (let i = 0; i < type_hash_list.length; i++) {
                        hash_list_len = hash_list.push(type_hash_list[i][1]);
                    }

                    let hash_query = {"hashes": hash_list}
                    let stmt_promise = getStatementByHash(hash_query)
                    
                    stmt_promise.then(function(stmt_response){
                        // Get statements array
                        var stmts = stmt_response.statements

                        // We could send an array of statement jsons, but then we 
                        // would have to keep track of which uuid is with which statement
                        // because I don't know if they're being returned in the same order
                        // as they were sent in. Instead, let's loop over statements and
                        // IDs for now

                        // Array to store query responses and corresponding uuids
                        var eng_res_array = [];
                        var stmt_uuid_array = [];

                        // Loop statements and store uuid and plain english query response
                        for (let i = 0; i < stmts.length; i++) {
                            stmt_json = stmts[i]
                            n = stmt_uuid_array.push(stmt_json.id)
                            json_stmt_array = {"statements": [stmt_json]}
                            eng_res_array.push(getEnglishByJson(json_stmt_array))
                        }

                        // Array Promises
                        // https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Promise/all
                        Promise.all(eng_res_array).then(function(eng_array) {

                            // Loop response array for plain english
                            for (let k = 0; k < eng_array.length; k++) {
                                sid = stmt_uuid_array[k]
                                var stmt = stmts[k]
                                eng_res = eng_array[k]
                                eng_plain = eng_res.sentences[sid]

                                // Output for Plain English
                                let AB_output_element_pe = document.createElement("h5")
                                AB_output_element_pe.textContent = "Statement " + k + ": " + eng_plain //+ "." // Add if plain english statement does not end in period
                                output_AB.appendChild(AB_output_element_pe)

                                // Loop evidence array
                                for (let j = 0; j < stmt.evidence.length; j++) {
                                    _pmid = stmt.evidence[j].pmid

                                    // Output for pmid link
                                    let AB_output_element_pmid = document.createElement("a")
                                    AB_output_element_pmid.href = "https://www.ncbi.nlm.nih.gov/pubmed/?term=" + _pmid
                                    AB_output_element_pmid.textContent = "PMID: " + _pmid
                                    output_AB.appendChild(AB_output_element_pmid)

                                    // Ouput for all evidence of of current stmt
                                    let AB_output_element_ev = document.createElement("div")
                                    AB_output_element_ev.textContent = "\"" + stmt.evidence[j].text.italics() + "\""
                                    output_AB.appendChild(AB_output_element_ev)
                                }
                            }
                        });
                    })
                },
                error: function() {
                    var output_AB = $("#expl_A_to_B")[0];
                    output_AB.innerHTML = null;
                    let AB_output_element_err = document.createElement("div")
                    AB_output_element_err.textContent = "Could not query at " + geneA_is_subj_address
                    output_AB.appendChild(AB_output_element_err)
                }

            })

            // Query and output all subj:B -> obj:A
            var geneB_is_subj_promise = $.ajax({
                url: geneB_is_subj_address,
                success: function(res) {
                    let obj = geneA
                    type_hash_list = res[obj]

                    var output_BA = $("#expl_B_to_A")[0];
                    output_BA.innerHTML = null;

                    // Create hash list
                    var hash_list = [];
                    var hash_list_len = 0;

                    for (let i = 0; i < type_hash_list.length; i++) {
                        hash_list_len = hash_list.push(type_hash_list[i][1]);
                    }

                    let hash_query = {"hashes": hash_list}
                    let stmt_promise = getStatementByHash(hash_query)

                    stmt_promise.then(function(stmt_response) {
                        var stmts = stmt_response.statements

                        var eng_res_array = [];
                        var stmt_uuid_array = [];

                        for (let i = 0; i < stmts.length; i++) {
                            stmt_json = stmts[i]
                            n = stmt_uuid_array.push(stmt_json.id)
                            json_stmt_array = {"statements": [stmt_json]}
                            eng_res_array.push(getEnglishByJson(json_stmt_array))
                        }

                        Promise.all(eng_res_array).then(function(eng_array) {

                            // Loop response array for plain english
                            for (let k = 0; k < eng_array.length; k++) {
                                sid = stmt_uuid_array[k]
                                var stmt = stmts[k]
                                eng_res = eng_array[k]
                                eng_plain = eng_res.sentences[sid]

                                // Output for Plain English
                                let BA_output_element_pe = document.createElement("h5")
                                BA_output_element_pe.textContent = "Statement " + k + ": " + eng_plain //+ "." // Add if plain english statement does not end in period
                                output_BA.appendChild(BA_output_element_pe)

                                // Loop evidence array
                                for (let j = 0; j < stmt.evidence.length; j++) {
                                    _pmid = stmt.evidence[j].pmid

                                    // Output for pmid link
                                    let BA_output_element_pmid = document.createElement("a")
                                    BA_output_element_pmid.href = "https://www.ncbi.nlm.nih.gov/pubmed/?term=" + _pmid
                                    BA_output_element_pmid.textContent = "PMID: " + _pmid
                                    output_BA.appendChild(BA_output_element_pmid)

                                    // Ouput for all evidence of of current stmt
                                    let BA_output_element_ev = document.createElement("div")
                                    BA_output_element_ev.textContent = stmt.evidence[j].text
                                    output_BA.appendChild(BA_output_element_ev)
                                }
                            }
                        })

                    })

                },
                error: function() {
                    var output_BA = $("#expl_B_to_A")[0];
                    output_BA.innerHTML = null;
                    let BA_output_element_err = document.createElement("div")
                    BA_output_element_err.textContent = "Could not query at " + geneB_is_subj_address
                    output_BA.appendChild(BA_output_element_err)
                }
            })

        } // Closing bracket for onChange: function()
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

            dropdownParent: "body",

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
                s3_prefix = "https://s3.amazonaws.com/depmap-public/prior_filtered_neighbor_lookup/neighbors_to_";
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

}); // This closes the load everything bracket
