// By putting the variable defintion here, it can be accessed easily in the console
var uuid_stmtjson_dict = {}; // used for buttons to be able to access resolved evidence etc

$(function(){
    // Load everything

    // Populate first dropdown from json array at:
    // https://s3.amazonaws.com/depmap-public/explainable_ids_1534216288.json
    var first_select_list = "https://s3.amazonaws.com/depmap-public/explainable_ids_1534216288.json";
    var select_first_gene, $select_first_gene
    var select_second_gene, $select_second_gene
    var indra_curation_addr = "https://lsm6zea7gg.execute-api.us-east-1.amazonaws.com/production/statements/from_hash/";
    // var indra_curation_addr = "http://127.0.0.1:5000/statements/from_hash/";
    var indra_server_addr = "https://lsm6zea7gg.execute-api.us-east-1.amazonaws.com/production/statements/from_hashes";
    // var indra_server_addr = "https://l3zhe2uu9c.execute-api.us-east-1.amazonaws.com/dev/statements/from_hashes";
    var indra_english_asmb = "http://api.indra.bio:8000/assemblers/english";

    // "GLOBAL" VARIABLE SCOPE
    var old_geneA = "A"
    var geneA = "A"
    var old_geneB = "B"
    var geneB = "B"
    // var uuid_stmtjson_dict = {}; // used for buttons to be able to access resolved evidence etc

    // OUTPUT POINTERS
    // Names
    // Complex
    var Aname_complex = document.getElementById("A_complex");
    var Bname_complex = document.getElementById("B_complex");
    // AB
    var Aname_AtoB = document.getElementById("A_AtoB");
    var Bname_AtoB = document.getElementById("B_AtoB");
    // BA
    var Aname_BtoA = document.getElementById("A_BtoA");
    var Bname_BtoA = document.getElementById("B_BtoA");
    // AXB
    var Aname_AXB = document.getElementById("A_AXB");
    var Bname_AXB = document.getElementById("B_AXB");
    // BXA
    var Aname_BXA = document.getElementById("A_BXA");
    var Bname_BXA = document.getElementById("B_BXA");
    // ABx
    var Aname_ABtoX = document.getElementById("A_ABtoX");
    var Bname_ABtoX = document.getElementById("B_ABtoX");
    // xAB
    var Aname_A_XtoAB = document.getElementById("A_XtoAB");
    var Bname_XtoAB = document.getElementById("B_XtoAB");

    // Output areas
    // Complex AB
    var output_AcB = $("#expl_A_complex_B")[0];
    // var output_ABcomplex = $("#AB_output_complex")[0];
    var AcB_ev_count = document.getElementById("collapseAcB_ev_count");
    // AB
    var output_AB = $("#expl_A_to_B")[0];
    // var output_AB_AB = $("#AB_output_AB")[0];
    var AB_ev_count = document.getElementById("collapseAB_ev_count");

    // BA
    var output_BA = $("#expl_B_to_A")[0];
    // var output_BA_BA = $("#BA_output_BA")[0];
    var BA_ev_count = document.getElementById("collapseBA_ev_count");

    // AXB
    var output_AXB = $("#expl_A_to_X_to_B")[0];
    var AXB_dd_div = "AXB_dropdown";
    var output_AX_AXB = $("#AX_output_AXB")[0];
    var output_XB_AXB = $("#XB_output_AXB")[0];
    var AXB_ev_count = document.getElementById("collapseAXB_ev_count");

    // BXA
    var output_BXA = $("#expl_B_to_X_to_A")[0];
    var BXA_dd_div = "BXA_dropdown";
    var output_BX_BXA = $("#BX_output_BXA")[0];
    var output_XB_BXA = $("#XA_output_BXA")[0];
    var BXA_ev_count = document.getElementById("collapseBXA_ev_count");

    // ABx
    var output_ABx = $('#expl_x_is_downstream')[0];
    var ABtox_dd_div = "ABtoX_dropdown";
    var output_AX_ABtoX = $("#AX_output_ABtoX")[0];
    var output_XB_ABtoX = $("#XB_output_ABtoX")[0];
    var ABx_ev_count = document.getElementById("collapse_st_X_count");

    // xAB
    var output_xAB = $('#expl_x_is_upstream')[0];
    var XtoAB_dd_div = "XtoAB_dropdown";
    var output_AX_XtoAB = $("#AX_output_XtoAB")[0];
    var output_XB_XtoAB = $("#XB_output_XtoAB")[0];
    var xAB_ev_count = document.getElementById("collapse_sr_X_count");

    function resetNamesOutput() {
        // console.log("Names reset");

        // Names for "A" and "B"
        Aname_complex.textContent = "A";
        Bname_complex.textContent = "B";
        Aname_AtoB.textContent = "A";
        Bname_AtoB.textContent = "B";
        Aname_BtoA.textContent = "A";
        Bname_BtoA.textContent = "B";
        Aname_AXB.textContent = "A";
        Bname_AXB.textContent = "B";
        Aname_BXA.textContent = "A";
        Bname_BXA.textContent = "B";
        Aname_ABtoX.textContent = "A";
        Bname_ABtoX.textContent = "B";
        Aname_A_XtoAB.textContent = "A";
        Bname_XtoAB.textContent = "B";

        // Clean up the output areas so old output doesn't stick around

        // Complex AB
        output_AcB.innerHTML = null;
        // output_ABcomplex.innerHTML = null;
        AcB_ev_count.textContent = "Statements: 0";
        AcB_ev_count.style = "background-color:#BBB; color: #FFFFFF;";
        // A->B
        output_AB.innerHTML = null;
        // output_AB_AB.innerHTML = null;
        AB_ev_count.textContent = "Statements: 0";
        AB_ev_count.style = "background-color:#BBB; color: #FFFFFF;";
        // B->A
        output_BA.innerHTML = null;
        // output_BA_BA.innerHTML = null;
        BA_ev_count.textContent = "Statements: 0";
        BA_ev_count.style = "background-color:#BBB; color: #FFFFFF;";
        // A->X->B
        output_AX_AXB.innerHTML = null;
        output_XB_AXB.innerHTML = null;
        AXB_ev_count.textContent = "X: 0";
        AXB_ev_count.style = "background-color:#BBB; color: #FFFFFF;";
        // B->X->A
        output_BX_BXA.innerHTML = null;
        output_XB_BXA.innerHTML = null;
        BXA_ev_count.textContent = "X: 0";
        BXA_ev_count.style = "background-color:#BBB; color: #FFFFFF;";
        // A<-X->B
        output_AX_XtoAB.innerHTML = null;
        output_XB_XtoAB.innerHTML = null;
        xAB_ev_count.textContent = "X: 0";
        xAB_ev_count.style = "background-color:#BBB; color: #FFFFFF;";
        // A->X<-B
        output_AX_ABtoX.innerHTML = null;
        output_XB_ABtoX.innerHTML = null;
        ABx_ev_count.textContent = "X: 0";
        ABx_ev_count.style = "background-color:#BBB; color: #FFFFFF;";
    }

    function allAreComplex(stmts) {
        for (hash of Object.keys(stmts)) {
            if (stmts[hash].type != "Complex") {
                return false;
            }
        }
        return true;
    }

    function sortByCol(arr, colIndex){
        arr.sort(sortFunction)
        function sortFunction(a, b) {
            a = a[colIndex]
            b = b[colIndex]
            return (a === b) ? 0 : (a < b) ? -1 : 1
        }
    }

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

    function getStatementsByHash(indra_query) {
        var api_key = document.getElementById("api_key_input").value;
        _url = indra_server_addr + "?api-key=" + api_key;
        stmts_db = $.ajax({
            url: _url,
            type: "POST",
            dataType: "json",
            contentType: "application/json",
            data: JSON.stringify(indra_query),
            });
        return stmts_db;
    };

    function getCurationHTMLByHash(hash) {
        ev_limit = 3;
        options = "&ev_limit=" + ev_limit
        url = indra_curation_addr + hash + "?format=html" + options
        return $.ajax({url: url})
    }

    function populateOutputDiv(outputDiv, curationHTMLtext) {
        fullHTML = document.createElement('html');
        fullHTML.innerHTML = curationHTMLtext;

        // Get div-class == "statement" and div-class == "evidence"
        statement = fullHTML.getElementsByClassName("statement")[0];
        evidence = fullHTML.getElementsByClassName("evidence")[0];

        // Append to output
        outputDiv.appendChild(statement);
        outputDiv.appendChild(evidence);
    }

    function grabJSON (url, callback) {
        return $.ajax({url: url, dataType: "json"});
    };

    // MOUSE HOVER LOAD PAGE
    $(".tiptext").mouseover(function() {
        $(this).children(".description").show();
    }).mouseout(function() {
        $(this).children(".description").hide();
    });

    $select_second_gene = $("#select_second_gene").selectize({
        // valueField: "second_item",
        valueField: "correlation",
        labelField: "second_item",
        searchField: ["second_item"],

        // A single field or an array of fields to sort by.
        sortField: {
            // field: "second_item",
            field: "correlation",
            direction: "desc"
        },

        onChange: function(value) {
            geneB = document.getElementById("select_second_gene").textContent.split(":")[0];

            // Reset Output and Names
            resetNamesOutput();

            // SET ADDRESSES TO AWS S3 DATA
            // Query of evidence for A->B
            // s3 bucket prefix string
            // Examples:
            // https://s3.amazonaws.com/depmap-public/indra_db_20180730_hash_json/CSDE1_is_subj.json <-- lookup objects given subject
            // https://s3.amazonaws.com/depmap-public/indra_db_20180730_hash_json/CSDE1_is_obj.json <-- lookup subjects given object
            // https://s3.amazonaws.com/depmap-public/correlation_pairs_above_03/correlates_with_A1BG.json
            // https://s3.amazonaws.com/depmap-public/Q3_depmap_20180730_db_explained_improved/A1BG_is_subj.json

            let s3_prefix = "https://s3.amazonaws.com/depmap-public/";
            // let s3_subj_expl = "Q3_depmap_20180730_db_explained_improved/";
            let s3_subj_expl = "pre_release_Q3_depmap_20180730_db_explained_improved/";
            let s3_indra_db = "indra_db_20180730_hash_json/"; // INDRA DB LOOKUP
            let s3_correlations = "correlation_pairs_above_03/correlates_with_";

            // Get the current selections
            geneA = document.getElementById("select_first_gene").value;
            geneB = document.getElementById("select_second_gene").textContent.split(":")[0];

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
                    //correlation_output.innerHTML = null;
                    //var correlation_output_element = document.createElement("a")
                    //var linkText = document.createTextNode("Link to depmap plot")
                    //correlation_output_element.appendChild(linkText);
                    //correlation_output_element.title = "Link to depmap plot for " + geneA + " vs " + geneB
                    //correlation_output_element.href = depmap1 + geneA + depmap2 + geneB + depmap3
                    correlation_output.href = depmap1 + geneA + depmap2 + geneB + depmap3
                    correlation_output.class = "nav-link active"
                    if (!isNumeric(correlation_AB)) {
                        // If this happens, you probably need to update the correlation jsons
                        console.log('Correlation is not a valid number!')
                        // correlation_output_element.textContent = "Link to depmap plot for " + geneA + " vs " + geneB + " not available."
                    }
                    // correlation_output.appendChild(correlation_output_element)
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
            var AcB_d_output = false;
            var AB_st_output = false;
            var AB_sr_output = false;
            var calledPMIDtitle = false;

            // Query and output all subj:A -> obj:B
            var geneA_is_subj_promise = $.ajax({
                url: geneA_is_subj_expl_address,
                success: function(res) {
                    let obj = geneB

                    // First to resolve: just set to true
                    if (!calledPMIDtitle) {
                        calledPMIDtitle = true;
                    } else {
                        var readyStateCheckInterval = setInterval(function() {
                            console.log("Current document readyState: " + document.readyState)

                            // options are "loading", "interactive", "complete"
                            if (document.readyState === "complete") {
                                clearInterval(readyStateCheckInterval);

                                // Source the curation row toggle function
                                $(function() {
                                    console.log("Re-reading toggle function")
                                    $("td[class='curation_toggle']").click(function(event) {
                                        console.log("Curation toggle click")
                                        event.stopPropagation();
                                        var $target = $(event.target);
                                        console.log($(event.target))
                                        if (event.target.dataset.clicked == "true") {
                                            console.log('trying to animate')
                                            // Toggle (animation duration in msec)
                                            $target.closest("tr").next().find("div").slideToggle(200);
                                        // First click event
                                        } else {
                                            console.log('first click event')
                                            // Stay down (animation duration in msec)
                                            $target.closest("tr").next().find("div").slideDown(400);

                                            // Change color of icon to light gray
                                            event.target.style="color:#A4A4A4;"

                                            // Set clicked to true
                                            event.target.dataset.clicked = "true"
                                        }
                                    });
                                });
                            }
                        }, 1000); // Time interval in milliseconds
                    }

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

                    // if connection undirected and not already printed
                    if (!AcB_d_output) {
                        if (connection_type_list.undirected.length > 0) {
                            var debug_string = 'output_AcB 1' // Kept for now in anticipation of future debugging needs

                            // Flag found so we don't make same call again for B->A
                            AcB_d_output = true

                            // Set names COMPLEX
                            Aname_complex.textContent = geneA;
                            Bname_complex.textContent = geneB;

                            // output_directs(output_pointer, source_output_pointer, ev_counter_pointer, type_hash_array, subj, obj, debug_string)
                            // output_directs(output_AcB, output_ABcomplex, AcB_ev_count, connection_type_list.undirected, geneA, geneB, debug_string);
                            output_directs(output_AcB, output_AcB, AcB_ev_count, connection_type_list.undirected, geneA, geneB, debug_string);
                        }
                    }

                    // if connection directed
                    if (connection_type_list.directed.length > 0) {
                        var debug_string = 'output_AB'

                        // Set names DIRECTed
                        Aname_AtoB.textContent = geneA;
                        Bname_AtoB.textContent = geneB;

                        // output_directs(output_pointer, source_output_pointer, ev_counter_pointer, type_hash_array, subj, obj, debug_string)
                        // output_directs(output_AB, output_AB_AB, AB_ev_count, connection_type_list.directed, geneA, geneB, debug_string);
                        output_directs(output_AB, output_AB, AB_ev_count, connection_type_list.directed, geneA, geneB, debug_string);
                    }

                    // 'x_is_intermediary'; This is for A->X->B
                    if (connection_type_list.x_is_intermediary.length > 0) {
                        var debug_string = 'output_AXB'

                        // Set names
                        Aname_AXB.textContent = geneA;
                        Bname_AXB.textContent = geneB;

                        // output_intermediary_new(output_pointer, SX_output_pointer, XO_output_pointer, x_counter_pointer, dd_select_id, x_array, geneA, geneB, geneA_lookup_address, geneB_lookup_address, debug_string)
                        output_intermediary_new(output_AXB, output_AX_AXB, output_XB_AXB, AXB_ev_count, AXB_dd_div, connection_type_list.x_is_intermediary, geneA, geneB, geneA_is_subj_address, geneB_is_obj_address, debug_string)
                    }

                    // 'x_is_downstream'
                    // Check if any output already
                    if (!AB_st_output) {
                        if (connection_type_list.x_is_downstream.length > 0) {
                            var debug_string = 'output_ABx 1'

                            // Set names
                            Aname_ABtoX.textContent = geneA;
                            Bname_ABtoX.textContent = geneB;

                            AB_st_output = true;

                            // output_intermediary_new(output_pointer, SX_output_pointer, XO_output_pointer, x_counter_pointer, dd_select_id, x_array, geneA, geneB, geneA_lookup_address, geneB_lookup_address, debug_string)
                            output_intermediary_new(output_ABx, output_AX_ABtoX, output_XB_ABtoX, ABx_ev_count, ABtox_dd_div, connection_type_list.x_is_downstream, geneA, geneB, geneA_is_subj_address, geneB_is_subj_address, debug_string)
                        }
                    }

                    // 'x_is_upstream':
                    // Check if any output already
                    if (!AB_sr_output) {
                        if (connection_type_list.x_is_upstream.length > 0) {
                            var debug_string = 'output_xAB 1'

                            // Flag found so we don't make same call again for B->A
                            AB_sr_output = true;

                            // Set names
                            Aname_A_XtoAB.textContent = geneA;
                            Bname_XtoAB.textContent = geneB;

                            // output_intermediary_new(output_pointer, SX_output_pointer, XO_output_pointer, x_counter_pointer, dd_select_id, x_array, geneA, geneB, geneA_lookup_address, geneB_lookup_address, debug_string)
                            output_intermediary_new(output_xAB, output_AX_XtoAB, output_XB_XtoAB, xAB_ev_count, XtoAB_dd_div, connection_type_list.x_is_upstream, geneA, geneB, geneA_is_obj_address, geneB_is_obj_address, debug_string)
                        }
                    }
                },
                error: function() {
                    output_AB.innerHTML = null;
                    let AB_output_element_err = document.createElement("div");
                    AB_output_element_err.textContent = "Could not query " + geneA_is_subj_expl_address;
                    output_AB.appendChild(AB_output_element_err);
                }

            })

            // Query and output all subj:B -> obj:A
            var geneB_is_subj_promise = $.ajax({
                url: geneB_is_subj_expl_address,
                success: function(res) {
                    let obj = geneA

                    // First to resolve: just set to true
                    if (!calledPMIDtitle) {
                        calledPMIDtitle = true;
                    } else {
                        var readyStateCheckInterval = setInterval(function() {
                            // options are "loading", "interactive", "complete"
                            if (document.readyState === "complete") {
                                clearInterval(readyStateCheckInterval);

                                // Source the curation row toggle function
                                $(function() {
                                    console.log("Re-reading toggle function")
                                    $("td[class='curation_toggle']").click(function(event) {
                                        console.log("Curation toggle click")
                                        event.stopPropagation();
                                        var $target = $(event.target);
                                        if (event.target.dataset.clicked == "true") {
                                            // Toggle (animation duration in msec)
                                            console.log('trying to animate')
                                            $target.closest("tr").next().find("div").slideToggle(200);
                                        // First click event
                                        } else {
                                            console.log('first click event')
                                            // Stay down (animation duration in msec)
                                            $target.closest("tr").next().find("div").slideDown(400);

                                            // Change color of icon to light gray
                                            event.target.style="color:#A4A4A4;"

                                            // Set clicked to true
                                            event.target.dataset.clicked = "true"
                                        }
                                    });
                                });
                            }
                        }, 1000); // Time interval in milliseconds
                    }

                    // Should return a dict of the format below
                    connection_type_list = res[obj]

                    // if connection undirected and not already printed
                    if (!AcB_d_output) {
                        if (connection_type_list.undirected.length > 0) {
                            AcB_d_output = true;
                            var debug_string = 'output_AcB 2'

                            // Set names COMPLEX
                            Aname_complex.textContent = geneA;
                            Bname_complex.textContent = geneB;

                            // output_directs(output_pointer, source_output_pointer, ev_counter_pointer, type_hash_array, subj, obj, debug_string)
                            // output_directs(output_AcB, output_ABcomplex, AcB_ev_count, connection_type_list.undirected, geneA, geneB, debug_string)
                            output_directs(output_AcB, output_AcB, AcB_ev_count, connection_type_list.undirected, geneA, geneB, debug_string)
                        }
                    }

                    // if connection directed
                    if (connection_type_list.directed.length > 0) {
                        var debug_string = 'output_BA'

                        // Set names DIRECTed
                        Aname_BtoA.textContent = geneA;
                        Bname_BtoA.textContent = geneB;

                        // output_directs(output_pointer, source_output_pointer, ev_counter_pointer, type_hash_array, subj, obj, debug_string)
                        // output_directs(output_BA, output_BA_BA, BA_ev_count, connection_type_list.directed, geneB, geneA, debug_string)
                        output_directs(output_BA, output_BA, BA_ev_count, connection_type_list.directed, geneB, geneA, debug_string)
                    }

                    // 'x_is_intermediary'; B->X->A
                    if (connection_type_list.x_is_intermediary.length > 0) {
                        var debug_string = 'output_BXA'

                        // Set names
                        Aname_BXA.textContent = geneA;
                        Bname_BXA.textContent = geneB;

                        // output_intermediary_new(output_pointer, SX_output_pointer, XO_output_pointer, x_counter_pointer, dd_select_id, x_array, geneA, geneB, geneA_lookup_address, geneB_lookup_address, debug_string)
                        output_intermediary_new(output_BXA, output_BX_BXA, output_XB_BXA, BXA_ev_count, BXA_dd_div, connection_type_list.x_is_intermediary, geneB, geneA, geneB_is_subj_address, geneA_is_obj_address, debug_string)
                    }

                    // Check if any output already
                    if (!AB_st_output) {
                        // 'x_is_downstream'
                        if (connection_type_list.x_is_downstream.length > 0) {
                            AB_st_output = true;
                            var debug_string = 'output_BAx 2'

                            // Set names
                            Aname_ABtoX.textContent = geneA;
                            Bname_ABtoX.textContent = geneB;

                            // output_intermediary_new(output_pointer, SX_output_pointer, XO_output_pointer, x_counter_pointer, dd_select_id, x_array, geneA, geneB, geneA_lookup_address, geneB_lookup_address, debug_string)
                            output_intermediary_new(output_ABx, output_AX_ABtoX, output_XB_ABtoX, ABx_ev_count, ABtox_dd_div, connection_type_list.x_is_downstream, geneA, geneB, geneA_is_subj_address, geneB_is_subj_address, debug_string)
                        }
                    }

                    // Check if any output already
                    if (!AB_sr_output) {
                        // 'x_is_upstream':
                        if (connection_type_list.x_is_upstream.length > 0) {
                            AB_sr_output = true;

                            var debug_string = 'output_xBA 2'

                            // Set names
                            Aname_A_XtoAB.textContent = geneA;
                            Bname_XtoAB.textContent = geneB;

                            // output_intermediary_new(output_pointer, SX_output_pointer, XO_output_pointer, x_counter_pointer, dd_select_id, x_array, geneA, geneB, geneA_lookup_address, geneB_lookup_address, debug_string)
                            output_intermediary_new(output_xAB, output_AX_XtoAB, output_XB_XtoAB, xAB_ev_count, XtoAB_dd_div, connection_type_list.x_is_upstream, geneA, geneB, geneA_is_obj_address, geneB_is_obj_address, debug_string)
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

            // An array of property names to analyze when filtering options.
            searchField: ["item"],

            // A single field or an array of fields to sort by.
            sortField: {
                field: "item",
                direction: "asc" 
            },

            // dropdownParent: "body",

            // Updates the current selection of first gene
            onChange: function(value) {
                if (!value.length) return;

                // Disable and clear second dropdown
                select_second_gene.clear();
                select_second_gene.clearOptions();
                select_second_gene.disable();

                geneB = "B"
                geneA = value

                // Reset output and names
                resetNamesOutput();

                // Set second query address example:
                // https://s3.amazonaws.com/depmap-public/prior_filtered_neighbor_lookup/neighbors_to_BRCA1.json
                // https://s3.amazonaws.com/depmap-public/neighbor_lookup/neighbors_to_A1BG.json
                // s3_prefix = "https://s3.amazonaws.com/depmap-public/neighbor_lookup/neighbors_to_"; // OLD
                s3_prefix = "https://s3.amazonaws.com/depmap-public/neighbor_lookup_2018sep19/neighbors_to_"; // NEW
                var second_dd_address = s3_prefix + value + ".json"

                // Query for next dropdown
                select_second_gene.load(function(callback) {
                    var second_json = $.ajax({
                        url: second_dd_address,
                        success: function(results) {
                            // var second_items = results.map(function(x) { return {second_item: x}; })
                            var second_items = results.map(function(x) { return {second_item: x[0] + ": correlation " + parseFloat(x[1]).toFixed(3).toString(), name: x[0], correlation: Math.abs(x[1]) }; })
                            select_second_gene.enable();
                            callback(second_items);
                        },
                        error: function() {
                            console.log("Could not load from " + second_dd_address)
                        }
                    })
                });
            } // This is closing bracket for "onChange: function(value)"
        }); // This is closing bracket for "$("#select_first_gene").selectize"

        var select_first_gene = $select_first_gene[0].selectize;
        var select_second_gene = $select_second_gene[0].selectize;

    });

    // Function for quering and outputting plain english description and statement evidence with PMIDs
    function output_directs(output_pointer, source_output_pointer, ev_counter_pointer, type_hash_array, subj, obj, debug_string){
        // console.log("< < Entering new output_directs call > >")

        // Create array to store each statement hash
        var hash_list = [];

        // type_hash_array contains [["type1", hash1], ["type2", hash2], ...]
        for (let i = 0; i < type_hash_array.length; i++) {
            hash_list.push(type_hash_array[i][1]);
        }
        nStmts = hash_list.length;
        ev_counter_pointer.textContent = "Statements: " + nStmts;
        if (nStmts == 0) {
            // No statements found, light gray
            ev_counter_pointer.style = "background-color:#BBB; color: #FFFFFF;"
        } else {
            // Statements found, set to standard gray 
            ev_counter_pointer.style = "background-color:#777; color: #FFFFFF;"
        }

        curationHTMLpromiseArray = [];
        for (let hash of hash_list) {
            curationHTMLpromiseArray.push(getCurationHTMLByHash(hash))
        }

        Promise.all(curationHTMLpromiseArray).then(function(curationHTMLtextArray) {
            for (let curationHTMLtext of curationHTMLtextArray) {
                populateOutputDiv(output_pointer, curationHTMLtext)
            }
        })
    } // Closes the output_directs function bracket

    // Use this function for A-X-B (same for all four) the query needs to be over two json lookups: SUBJ_is_subj and OBJ_is_obj
    function output_intermediary_new(output_pointer, SX_output_pointer, XO_output_pointer, x_counter_pointer, dd_select_id, x_array, geneA, geneB, geneA_lookup_address, geneB_lookup_address, debug_string){
        // console.log(('Called output_intermediary_new from ' + debug_string))
        var items = x_array.map(function(x) { return { x_value: x[0], item: x[0] + ": rank " + parseFloat(x[1]).toFixed(3).toString(), rank: x[1] }; })

        // Update the count of X in the badge
        x_counter_pointer.textContent = "X: " + x_array.length.toString()
        if (x_array.length == 0) {
            // No X found, gray out 
            x_counter_pointer.style = "background-color:#BBB; color: #FFFFFF;"
        } else {
            // X found, set to black
            x_counter_pointer.style = "background-color:#777; color: #FFFFFF;"
        }

        // Load list of X into dropdown; ALTERNATIVELY, MAKE DROPDOWN FUNCTIONS FOR EACH A-X-B EXPLICITLY
        $select_intermediate = $("#"+dd_select_id).selectize({
            options: items,
            valueField: "x_value",
            labelField: "item",
            searchField: ["item"],

            // A single field or an array of fields to sort by.
            sortField: {
                field: "rank",
                direction: "desc"
            },

            // General settings
            create: false, // Don't allow user to add options
            maxItems: 1, // Only allow one item to be selected
            placeholder: "Select gene X...", // Placeholder when nothing is selected

            // On select/change: Query A-X and B-X and output the english statements and their evidence
            // Also clear the output area so that a fresh one can be sent once a new X is selected
            onChange: function(x_value) {
                if (!x_value.length) return;
                two_promises = [];
                two_promises.push(grabJSON(geneA_lookup_address)) // <--- Query for 'SUBJ_is_subj'
                two_promises.push(grabJSON(geneB_lookup_address)) // <--- Query for 'OBJ_is_obj'

                // Wait for both promises to resolve
                Promise.all(two_promises).then(function(two_jsons_ar){
                    // Get the the hash arrays
                    geneA_lookup = two_jsons_ar[0];
                    geneB_lookup = two_jsons_ar[1];

                    // Create pointers to nothing so that we can give the x_counter_pointer something
                    let SX_fake_x_counter = document.createElement("div")
                    let XO_fake_x_counter = document.createElement("div")

                    // null the output area and create new <div>s for both outputs
                    // A-X OUTPUT
                    SX_output_pointer.innerHTML = null;
                    let SX_output_div = document.createElement("div")
                    let SX_output_header = document.createElement("h4")
                    SX_output_header.style = "background-color:#F2F2F2;"
                    SX_output_header.textContent = geneA + ", " + x_value;
                    SX_output_div.appendChild(SX_output_header)
                    SX_output_pointer.appendChild(SX_output_div)
                    output_pointer.appendChild(SX_output_pointer)
                    // X-B OUTPUT
                    XO_output_pointer.innerHTML = null;
                    let XO_output_div = document.createElement("div")
                    let XO_output_header = document.createElement("h4")
                    XO_output_header.style = "background-color:#F2F2F2;"
                    XO_output_header.textContent = x_value + ", " + geneB;
                    XO_output_div.appendChild(XO_output_header)
                    XO_output_pointer.appendChild(XO_output_div)
                    output_pointer.appendChild(XO_output_pointer)

                    // output_directs(output_pointer, source_output_pointer, ev_counter_pointer, type_hash_array, subj, obj, debug_string)
                    output_directs(SX_output_pointer, SX_output_div, SX_fake_x_counter, geneA_lookup[x_value], geneA, x_value, debug_string)
                    output_directs(XO_output_pointer, XO_output_div, XO_fake_x_counter, geneB_lookup[x_value], x_value, geneB, debug_string)
                });
            }
        })
    } // Closes the output_intermediary_new function bracket
}); // This closes the load everything bracket
