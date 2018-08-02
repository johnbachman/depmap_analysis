// Populate list from json list at:
// https://s3.amazonaws.com/depmap-public/all_hgnc_prior_ll03.json
all_hgnc_prior_ll03 = 'https://s3.amazonaws.com/depmap-public/all_hgnc_prior_ll03.json';

function grabJSON (url, callback) {
    return $.ajax({url: url,});
};

// Save response in a promise
var json_promise = grabJSON(all_hgnc_prior_ll03);

// Use the .then() of promise to make the call only after we have a response. By some JS magic, res == json_promise.responseJSON
json_promise.then(function(res) {
    data = res
    var items = data.map(function(x) { return { item: x }; })

    $('#select_subj').selectize({
        // USage documentation: https://github.com/selectize/selectize.js/blob/master/docs/usage.md

        // Populate the dropdown options from an array ('items')
        options: items,

        // The name of the property in 'items' to populate the dropdown with.
        labelField: 'item',

        // The name of the property to use as the value when an item is selected.
        valueField: 'item',

        // Allows the user to create new items not in the initial list.
        // create: true,

        // An array of property names to analyze when filtering options.
        searchField: ['item'],

        // A single field or an array of fields to sort by.
        sortField: {
            field: 'item',
            direction: 'asc' 
        },

        // Updates the current selection 
        onChange:function(value) {

            // Get the dropdown 'select' block
            let subj_input = document.getElementById('select_subj');
            
            // Refer to the div (or other object) where the output text should be
            var output_text = $('#my_output')[0];

            // Add an empty innerHTML object (otherwise it keeps appending to current HTML object)
            output_text.innerHTML = null;
            
            // build an element (here: 'thingy') in the innerHTML that includes the selected value.
            var thingy = document.createElement("span");

            // Get the text from the selected item dropdown 
            thingy.textContent = 'Subject: ' + subj_input.options[subj_input.selectedIndex].text

            // Append the element to the div object
            output_text.appendChild(thingy)
        },
        dropdownParent: 'body'
    });
});
