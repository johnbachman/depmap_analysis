const SUBMIT_URL = './query/submit';
const S3_QUERY_CAHCE = 'https://bigmech.s3.amazonaws.com/indra_network_search/';
const MAX_K_PATHS = 50;
let pathStmtHashes = [];

let stmtItems = [];
for (s of stmtOptions) {
  stmtItems.push({
    value: s.toLowerCase(),
    label: s,
    selected: false,
    disabled: false
  })
}

let nodeItems = [];
for (let n of nodeOptions) {
  let sel = false;
  if (n === 'HGNC') {
    sel = true;
  }
  nodeItems.push({
    value: n.toLowerCase(),
    label: n,
    selected: sel,
    disabled: false
  });
}

let termNamespaces = [];
for (let n of termNsOptions) {
  let sel = false;
  termNamespaces.push({
    value: n.toLowerCase(),
    label: n,
    selected: sel,
    disabled: false
  });
}

function submitQuery() {
  let beliefEntry = parseInRange(document.getElementById('belief-cutoff').value,
                                 document.getElementById('belief-cutoff').min,
                                 document.getElementById('belief-cutoff').max,
                                 false);
  let kShortestEntry = parseInRange(document.getElementById('k-shortest').value,
                                    document.getElementById('k-shortest').min,
                                    document.getElementById('k-shortest').max,
                                    true);
  let maxPerNode = parseInRange(document.getElementById('max-per-node').value,
                                document.getElementById('max-per-node').min,
                                document.getElementById('max-per-node').max,
                                true)
  let timeoutEntry = parseInRange(document.getElementById('user_timeout').value,
                                  document.getElementById('user_timeout').min,
                                  document.getElementById('user_timeout').max,
                                  false);
  let nodeFilterList = [];
  for (c of document.getElementById('node-filter').children) {
    nodeFilterList.push(c.value)
  }
  if (nodeFilterList.length === 0) {
    alert('Must select at least one node namespace to include in path');
    return;
  }
  let stmtFilterList = [];
  for (c of document.getElementById('stmt-filter').children) {
    stmtFilterList.push(c.value)
  }
  let termNsList = [];
  for (c of document.getElementById('terminal-namespaces').children) {
    termNsList.push(c.value)
  }
  if (!document.getElementById('fplx-edges').checked) {
    stmtFilterList.push('fplx')
  }
  let nodeBlackList = [];
  for (nn of document.getElementById('node-blacklist').value.split(',')) {
    // Strip leading and trailing whitespace and push to array
    if (nn.replace(/\s/g, '')) nodeBlackList.push(nn.replace(/\s/g, ''));
  }
  let edgeHashBlacklist = [];
  for (eh of document.getElementById('edge-hash-blacklist').value.split(',')) {
    // Strip whitespace
    if (eh.replace(/\s/g, '')) edgeHashBlacklist.push(eh.replace(/\s/g, ''))
  }
  let cullBestNode = parseInt(document.getElementById('cull-best-node').value) || 0;
  if (cullBestNode && cullBestNode < 1) {
    alert('Culling best node every # paths must be positive integer');
    return;
  }
  let meshIdList = []
  for (id of document.getElementById('mesh-id-list').value.split(',')) {
    // Strip whitespace
    if (id.replace(/\s/g, ''))  meshIdList.push(id.replace(/\s/g, ''))
  }
  let constC = parseInRange(document.getElementById('const_c').value,
                            document.getElementById('const_c').min,
                            document.getElementById('const_c').max,
                            true);
  let constTk = parseInRange(document.getElementById('const_tk').value,
                             document.getElementById('const_tk').min,
                             document.getElementById('const_tk').max,
                             true);
  let statusBox = document.getElementById('query-status');
  let source = document.getElementById('source').value;
  let target = document.getElementById('target').value;
  let queryDict = { 
    source: source,
    target: target,
    stmt_filter: stmtFilterList,
    edge_hash_blacklist: edgeHashBlacklist,
    node_filter: nodeFilterList,
    node_blacklist: nodeBlackList,
    path_length: parseInt(document.getElementById('path-length').value) || false,
    sign: document.getElementById('sign-dd').value,
    weighted: document.getElementById('weighted-search').checked,
    bsco: beliefEntry,
    // direct_only: document.getElementById('direct-only').checked,
    curated_db_only: document.getElementById('curated-db-only').checked,
    fplx_expand: document.getElementById('fplx-expand').checked,
    k_shortest: kShortestEntry,
    max_per_node: maxPerNode,
    cull_best_node: cullBestNode,
    mesh_ids: meshIdList,
    strict_mesh_id_filtering: document.getElementById('strict-mesh-id-filtering').checked,
    const_c: constC,
    const_tk: constTk,
    user_timeout: timeoutEntry,
    two_way: document.getElementById('two-ways').checked,
    shared_regulators: document.getElementById('shared-regulators').checked,
    terminal_ns: termNsList,
    format: 'html'
  };

  // Hide on submit
  document.getElementById('download-link').hidden = true;

  // Empty hash list on submit
  pathStmtHashes = [];

  let _url = SUBMIT_URL;
  statusBox.textContent = 'Query submitted...';
  console.log('Query submitted:');
  console.log(queryDict);
  let response = $.ajax({
    url: _url,
    type: 'POST',
    dataType: 'json',
    contentType: 'application/json',
    data: JSON.stringify(queryDict),
    complete: function(xhr) {
      console.log(xhr);
      switch (xhr.status) {
        case 200:
          if (xhr.responseText) {
            window.location.replace(xhr.responseText)
          }
          break;
        case 202:
          if (xhr.responseText) {
            console.log(xhr.responseText)
          }
          statusBox.textContent = xhr.responseText;
          break;
        case 504:
          // Server timeout
          statusBox.textContent = 'Error: 504: Gateway Time-out';
          // ToDo Add more messages for different HTML errors
          break;
        default:
          console.log('Submission error: check ajax response');
          statusBox.textContent = 'Error: ' + xhr.status + ': ' + xhr.responseText;
          break;
      }
    }

  });
  console.log(response);
}

function isEmptyResult(resultJson, allowTimeOut=false) {
  /* Check if the resulting json from the query is empty
  EMPTY_RESULT = {'paths_by_node_count': {'forward': {}, 'backward': {}},
                         'common_targets': [],
                         'common_parents': {},
                         'timeout': false};
  With 'allowTimeOut' set to true, only the result itself is considered,
  regardless of time out status. If 'allowTimeOut' is false, the empty result
  has to be from a non-timed out query to be considered empty.
  */
  if (allowTimeOut) {
  return $.isEmptyObject(resultJson.paths_by_node_count) &&
    $.isEmptyObject(resultJson.paths_by_node_count.forward) &&
    $.isEmptyObject(resultJson.paths_by_node_count.backward) &&
    $.isEmptyObject(resultJson.common_targets) &&
    $.isEmptyObject(resultJson.common_parents);
  } else {
  return $.isEmptyObject(resultJson.paths_by_node_count) &&
    $.isEmptyObject(resultJson.paths_by_node_count.forward) &&
    $.isEmptyObject(resultJson.paths_by_node_count.backward) &&
    $.isEmptyObject(resultJson.common_targets) &&
    $.isEmptyObject(resultJson.common_parents) &&
    !(resultJson.timeout);
  }
}

function fillOldQuery(oldQueryJson) {
  // Input boxes: mandatory
  document.getElementById('source').value = oldQueryJson.source;
  document.getElementById('target').value = oldQueryJson.target;
  // Input boxes: optional
  if (!($.isEmptyObject(oldQueryJson.edge_hash_blacklist))) {
    document.getElementById('edge-hash-blacklist').value = oldQueryJson.edge_hash_blacklist.join(', ')
  }
  if (!($.isEmptyObject(oldQueryJson.node_blacklist))) {
    document.getElementById('node-blacklist').value = oldQueryJson.node_blacklist.join(', ')
  }
  if (oldQueryJson.path_length) document.getElementById('path-length').value = parseInt(oldQueryJson.path_length);
  if (oldQueryJson.sign) document.getElementById('sign-dd').value = oldQueryJson.sign;
  if (oldQueryJson.bsco) document.getElementById('belief-cutoff').value = parseInRange(oldQueryJson.bsco, 0, 1, false);
  if (oldQueryJson.k_shortest) document.getElementById('k-shortest').value = parseInt(oldQueryJson.k_shortest);
  if (oldQueryJson.max_per_node) document.getElementById('max-per-node').value = parseInt(oldQueryJson.max_per_node);
  if (oldQueryJson.cull_best_node) {
    document.getElementById('cull-best-node').value = parseInt(oldQueryJson.cull_best_node)
  }
  if (oldQueryJson.user_timeout) {
    document.getElementById('user_timeout').value = parseInRange(oldQueryJson.user_timeout, 0, 999, false)
  }
  if (oldQueryJson.mesh_ids) document.getElementById('mesh-id-list').value = oldQueryJson.mesh_ids.join(', ');
  // Checkboxes
  if (oldQueryJson.shared_regulators) document.getElementById('shared-regulators').checked = oldQueryJson.shared_regulators;
  if (oldQueryJson.strict_mesh_id_filtering) {
    document.getElementById('strict-mesh-id-filtering').checked = oldQueryJson.strict_mesh_id_filtering
  }
  document.getElementById('const_c').value = parseInRange(oldQueryJson.const_c, 1, 120, true);
  document.getElementById('const_tk').value = parseInRange(oldQueryJson.const_tk, 1, 120, true);

  if (oldQueryJson.weighted) document.getElementById('weighted-search').checked = oldQueryJson.weighted;
  // if (oldQueryJson.direct_only) document.getElementById('direct-only').checked = oldQueryJson.direct_only;
  if (oldQueryJson.curated_db_only) document.getElementById('curated-db-only').checked = oldQueryJson.curated_db_only;
  if (oldQueryJson.fplx_expand) document.getElementById('fplx-expand').checked = oldQueryJson.fplx_expand;
  if (oldQueryJson.two_way) document.getElementById('two-ways').checked = oldQueryJson.two_way;
  if (!($.isEmptyObject(oldQueryJson.stmt_filter))) {
    document.getElementById('fplx-edges').checked = !oldQueryJson.stmt_filter.includes('fplx')
  }

  let stmtItems = [];
  let selStmts = oldQueryJson.stmt_filter;
  for (s of stmtOptions) {
    let sel = selStmts.includes(s.toLowerCase());
    stmtItems.push({
      value: s.toLowerCase(),
      label: s,
      selected: sel,
      disabled: false
    })
  }

  let nodeItems = [];
  let selNodes = oldQueryJson.node_filter;
  for (let n of nodeOptions) {
    let sel = selNodes.includes(n.toLowerCase());
    nodeItems.push({
      value: n.toLowerCase(),
      label: n,
      selected: sel,
      disabled: false
    })
  }

  let termNamespaces = [];
  let selTermNs = oldQueryJson.terminal_ns;
  for (let n of termNsOptions) {
    let sel = selTermNs.includes(n.toLowerCase());
    termNamespaces.push({
      value: n.toLowerCase(),
      label: n,
      selected: sel,
      disabled: false
    })
  }
  return [stmtItems, nodeItems, termNamespaces]
}

function fillResultsTable(data, source, target){
  console.log(data);
  const statusBox = document.getElementById('query-status');
  if (isEmptyResult(data.result, false) && !source && !target) {
    console.log('Empty result!')
    return false;
  }

  let downloadURL = `${S3_QUERY_CAHCE}${data.query_hash}_result.json`;
  let downloadLink = '';
  if (downloadURL) {
    downloadLink = ` Click <a href="${downloadURL}" download>here</a> to download the results as a json`
  }
  let idCounter = 0;
  let currentURL = window.location.href.split('#')[0]
  if ((data.result.common_targets && data.result.common_targets.length > 0) ||
      (data.result.common_parents && Object.keys(data.result.common_parents).length > 0) ||
      (data.result.paths_by_node_count && Object.keys(data.result.paths_by_node_count).length > 0)) {
    if (data.result.timeout) statusBox.innerHTML = 'Query timed out with' +
      ' results!' + downloadLink;
    else if (data.result.node_not_found) statusBox.innerHTML = 'Error: ' + data.result.node_not_found;
    else statusBox.innerHTML = 'Query resolved!' + downloadLink;
    // for each path length:
    //   for each path:
    //     Output: path | list supporting statements per edge
    let commonParents = data.result.common_parents;
    let pathsKeyedArrayForward = data.result.paths_by_node_count.forward;
    let pathsKeyedArrayBackward = data.result.paths_by_node_count.backward;
    let simpleCommonTargets = data.result.common_targets;
    let simpleSharedRegulators = data.result.shared_regulators;
    let tableArea = document.getElementById('table-area');
    pathStmtHashes = data.result.paths_by_node_count.path_hashes;

    // Fill common parents table
    if (commonParents &&
        commonParents.common_parents &&
        commonParents.common_parents.length > 0) {
      let cardHtml = generateCommonParents(currentURL);
      idCounter += 1;
      tableArea.appendChild(cardHtml);
      document.getElementById('subject-placeholder-cp').textContent =
        `${commonParents.source_id}@${commonParents.source_ns}`;
      document.getElementById('object-placeholder-cp').textContent =
        `${commonParents.target_id}@${commonParents.target_ns}`;
      let cpTableBody = document.getElementById('query-results-cp');
      document.getElementById('npaths-cp').textContent = 'Entities: ' + commonParents.common_parents.length;
      for (let par of commonParents.common_parents) {
        let newRow = document.createElement('tr');

        let newIdCol = document.createElement('td');

        let [parNS, parID, parUri] = par;
        let hs = `<a href="${parUri}" target="_blank">${parID}</a>`;
        newIdCol.innerHTML = hs;
        newRow.appendChild(newIdCol);

        let newNsCol = document.createElement('td');
        newNsCol.textContent = parNS.toUpperCase();
        newRow.appendChild(newNsCol);

        cpTableBody.appendChild(newRow)
      }
    }

    // Fill common targets table
    if (simpleCommonTargets && simpleCommonTargets.length > 0) {
      let cardHtml = generateCommonTargets();
      tableArea.appendChild(cardHtml);
      let ctTableBody = document.getElementById('query-results-common-targets');
      document.getElementById('common-targets').textContent = `Targets: ${simpleCommonTargets.length}`;
      document.getElementById('subject-placeholder-ct').textContent = source;
      document.getElementById('object-placeholder-ct').textContent = target;
      for (let targetDict of simpleCommonTargets) {
        for (let key in targetDict) {
          if (key !== 'lowest_highest_belief') {
            let newRow = document.createElement('tr');

            let newTargetCol = document.createElement('td');
            let linkId = `linkout${idCounter}`;
            newTargetCol.id = linkId;
            idCounter += 1;
            newTargetCol.innerHTML = `<a class="stmt_toggle" title="Right click and copy link to link to this result" href="${currentURL}#${linkId}">${key}</a>`;
            newRow.appendChild(newTargetCol);

            let newTargetPaths = document.createElement('td');
            newTargetPaths.innerHTML = generateTargetLinkout(targetDict[key]);
            newRow.appendChild(newTargetPaths);

            ctTableBody.appendChild(newRow);
          }
        }
      }
    }

    // Fill shared regulators
    if (simpleSharedRegulators && simpleSharedRegulators.length > 0){
      let cardHtml = generateSharedRegulators();
      tableArea.appendChild(cardHtml);
      let srTableBody = document.getElementById('query-results-shared-regulators');
      document.getElementById('shared-regulators-span').innerText = `Regulators: ${simpleSharedRegulators.length}`;
      document.getElementById('subject-placeholder-sr').textContent = source;
      document.getElementById('object-placeholder-sr').textContent = target;
      for (let regulatorDict of simpleSharedRegulators) {
        for (let key in regulatorDict) {
          if (key !== 'lowest_highest_belief') {
            let newRow = document.createElement('tr');

            let newRegulatorCol = document.createElement('td');
            let linkId = `linkout${idCounter}`;
            idCounter += 1;
            newRegulatorCol.id = linkId;
            newRegulatorCol.innerHTML = `<a class="stmt_toggle" title="Right click and copy link to link to this result" href="${currentURL}#${linkId}">${key}</a>`;
            newRow.appendChild(newRegulatorCol);

            let newRegulatorPaths = document.createElement('td');
            newRegulatorPaths.innerHTML = generateTargetLinkout(regulatorDict[key]);
            newRow.appendChild(newRegulatorPaths);

            srTableBody.appendChild(newRow)
          }
        }
      }

    };

    // Fill directed paths tables
    if (pathsKeyedArrayForward && Object.keys(pathsKeyedArrayForward).length > 0) {
      for (let len in pathsKeyedArrayForward) {
        let cardHtml = generateCardTable(len-1, 'f');
        tableArea.appendChild(cardHtml);
        // cardHtml.getElementById('npaths-' + (len-1)).textContent = 'Paths: ' + pathsKeyedArrayForward[len].length;
        cardHtml.getElementsByClassName('badge')[0].textContent = 'Paths: ' + pathsKeyedArrayForward[len].length;
        // cardHtml.getElementById('subject-placeholder-' + (len-1)).textContent = source;
        cardHtml.getElementsByClassName('subject-placeholder')[0].textContent = source;
        // cardHtml.getElementById('object-placeholder-' + (len-1)).textContent = target;
        cardHtml.getElementsByClassName('object-placeholder')[0].textContent = target;
        // let tableBody = cardHtml.getElementById('query-results-' + (len-1));
        let tableBody = cardHtml.getElementsByClassName('table-body')[0];
        for (let pathArray of pathsKeyedArrayForward[len]) {
          let newRow = document.createElement('tr');
          // Add current path
          let newPath = document.createElement('td');
          let linkId = `linkout${idCounter}`;
          newPath.id = linkId;
          idCounter += 1;
          newPath.innerHTML = `<a class="stmt_toggle" title="Right click and copy link to link to this result" href="${currentURL}#${linkId}">${pathArray.path.join('&rarr;')}</a>`;
          newRow.appendChild(newPath);

          let newWeights = document.createElement('td');
          newWeights.innerHTML = generatePathWeights(pathArray.stmts);
          newRow.appendChild(newWeights)

          let newSupport = document.createElement('td');
          newSupport.innerHTML = generatePathLinkout(pathArray.stmts);
          newRow.appendChild(newSupport);

          tableBody.appendChild(newRow);
        }
      }
    }
    if (pathsKeyedArrayBackward && Object.keys(pathsKeyedArrayBackward).length > 0) {
      let tableArea = document.getElementById('table-area');
      for (len in pathsKeyedArrayBackward) {
        let cardHtml = generateCardTable(len-1, 'b');
        tableArea.appendChild(cardHtml);
        // cardHtml.getElementById('npaths-' + (len-1)).textContent = 'Paths: ' + pathsKeyedArrayBackward[len].length;
        cardHtml.getElementsByClassName('badge')[0].textContent = 'Paths: ' + pathsKeyedArrayBackward[len].length;
        // cardHtml.getElementById('subject-placeholder-' + (len-1)).textContent = target;
        cardHtml.getElementsByClassName('subject-placeholder')[0].textContent = target;
        // cardHtml.getElementById('object-placeholder-' + (len-1)).textContent = source;
        cardHtml.getElementsByClassName('object-placeholder')[0].textContent = source;
        // let tableBody = cardHtml.getElementById('query-results-' + (len-1));
        let tableBody = cardHtml.getElementsByClassName('table-body')[0];
        for (pathArray of pathsKeyedArrayBackward[len]) {
          let newRow = document.createElement('tr');
          // Add current path
          let newPath = document.createElement('td');
          let linkId = `linkout${idCounter}`;
          newPath.id = linkId;
          idCounter += 1;
          newPath.innerHTML = `<a class="stmt_toggle" title="Right click and copy link to link to this result" href="${currentURL}#${linkId}">${pathArray.path.join('&rarr;')}</a>`;
          newRow.appendChild(newPath);

          let newWeights = document.createElement('td');
          newWeights.innerHTML = generatePathWeights(pathArray.stmts);
          newRow.appendChild(newWeights)

          let newSupport = document.createElement('td');
          newSupport.innerHTML = generatePathLinkout(pathArray.stmts);
          newRow.appendChild(newSupport);

          tableBody.appendChild(newRow);
        }
      }
    }
  } else {
    if (data.result.timeout) statusBox.textContent = 'Query timed out!';
    else statusBox.textContent = 'No path found, try expanding the search parameters';
  }

  // If we have hashes, show download link
  if (pathStmtHashes.length > 0) {
    document.getElementById('download-link').href = `/stmts_download/stmts.json?query=${data.query_hash}`;
    document.getElementById('download-p').hidden = false;
  }
}

function clearAllTables() {
  // Delete cards
  $('.card').remove();
}

function resetCounters() {
  for (c of document.getElementsByClassName('path-count')) {
    let text = c.textContent.split(': ')[0];
    c.textContent = text + ': 0';
  }
}

function isNumeric(n) {
  return !isNaN(parseFloat(n)) && isFinite(n);
}

function parseInRange(num, min, max, int) {
  if (!num) return false;

  num = parseFloat(num);
  if (!num) return false;

  min = parseFloat(min);
  max = parseFloat(max);

  let ret = num;

  if (num < min) ret = min;
  else if (num > max) ret = max;

  if (int) ret = parseInt(ret);

  return ret;
}

// Function to put stmts from path in a HTML format with linkouts
function generatePathLinkout(pathArray) {
  // Build an html string of the paths
  // pathArray is a list of dicts for each edge in the path:
  // pathArray[i] gives a dict for the support of edge i
  // pathArray[i][stmt type j] is a list of all statements of type j supporting edge i
  // pathArray[i][stmt type j][k] is statement k of type j supporting edge i
  let htmlString = '';
  for (let edgeDict of pathArray) {
    if (Object.keys(edgeDict).length > 0) {
      let subj = edgeDict.subj;
      let obj = edgeDict.obj;
      htmlString += `<h5><b>${subj} &rarr; ${obj}</b><span class="float-right">Statement types: ${(Object.keys(edgeDict).length-2).toString()}</span></h5>`;
      for (let stmt_type in edgeDict) {
        if ((stmt_type !== 'subj') && (stmt_type !== 'obj') && (stmt_type !== 'weight_to_show')) {
          // let dbLink = '';
          // if (stmt.stmt_hash.startsWith('http')) dbLink = stmt.stmt_hash;
          // ?subject=BRCA1&object=BRCA2&type=activation&ev_limit=1
          let agentsString = '';
          if (stmt_type.toLowerCase() === 'complex') agentsString = `agent0=${subj}&agent1=${obj}&`;
          else agentsString = `subject=${subj}&object=${obj}&`;
          let dbLink = INDRA_DB_URL_AGENTS + agentsString + `type=${stmt_type}&ev_limit=1`;
          let sourceBadges = generateSourceBadges(edgeDict[stmt_type]);
          htmlString += '<a href="' + dbLink + '" target="_blank">' + subj + ', ' +
            stmt_type + ', ' + obj + '</a>' + sourceBadges + '<br>'
        }
      }
    }
  }
  return htmlString.substring(0, htmlString.length-4); // Cut out the last <br>
}

function generatePathWeights(pathArray) {
  let htmlString = '';
  for (let edgeDict of pathArray) {
    if (Object.keys(edgeDict).length > 0) {
      let subj = edgeDict.subj;
      let obj = edgeDict.obj;
      htmlString += '<h5>' + edgeDict.weight_to_show + '</h5>';
      for (let stmt_type in edgeDict)
        if ((stmt_type !== 'subj') && (stmt_type !== 'obj') && (stmt_type !== 'weight_to_show'))
          htmlString += '<a/><br>';
    }
  }
  return htmlString.substring(0, htmlString.length-4); // Cut out the last <br>
}

function generateSourceBadges(stmt_list) {
  // ToDo Add upp source counts in source dict
  let all_source_counts = {};
  for (let stmt of stmt_list) {
    let source_counts = stmt.source_counts;
    for (let source in source_counts) {
      if (all_source_counts[source]) all_source_counts[source] += source_counts[source];
      else all_source_counts[source] = source_counts[source];
    }
  }
  let htmlString = '';
  let badgeStart = '<span class="badge badge-secondary badge-pill float-right ';
  let badgeEnd = '</span>';
  for (let source in all_source_counts) {
    htmlString += badgeStart + 'source-'+ source + '">' + source + ': ' + all_source_counts[source].toString() + badgeEnd
  }
  return htmlString
}

function generateTargetLinkout(pathPair) {
  let htmlString = '';
  for (let paths of pathPair) {
    for (let path of paths) {
      let subj = path.subj;
      let obj = path.obj;
      htmlString += `<h5><b>${subj} &rarr; ${obj}</b><span class="float-right">Support: ${(Object.keys(path).length-2).toString()}</span></h5>`;
      for (let stmt_type in path) {
        if ((stmt_type !== 'subj') && (stmt_type !== 'obj') && (stmt_type !== 'weight_to_show')) {
          // let dbLink = '';
          // if (stmt.stmt_hash.startsWith('http')) dbLink = stmt.stmt_hash;
          // else dbLink = INDRA_DB_URL_HASH + stmt.stmt_hash + '?format=html';
          let agentsString = '';
          if (stmt_type.toLowerCase() === 'complex') agentsString = `agent0=${subj}&agent1=${obj}&`;
          else agentsString = `subject=${subj}&object=${obj}&`;
          let dbLink = INDRA_DB_URL_AGENTS + agentsString + `type=${stmt_type}&ev_limit=1`;
          let sourceBadges = generateSourceBadges(path[stmt_type]);
          htmlString += '<a href="' + dbLink + '" target="_blank">' + subj + ', ' +
            stmt_type + ', ' + obj + '</a>' + sourceBadges + '<br>'
        }
      }
    }
  }
  return htmlString.substring(0, htmlString.length-4); // Cut out the last <br>
}

function generateCardTable(len, dir) {
  let intermArrows = '';
  if (len > 1) {
    for (let i = 1; i < len; i++) {
      intermArrows += 'X' + i + '&rarr;'
    }
  }

  let newCard = document.createElement('div');
  newCard.className = 'card';
  newCard.innerHTML = '<div class="card-header" data-toggle="collapse" ' +
    'data-target="#collapse-paths-' + len + dir + ' " aria-expanded="true"' +
    ' aria-controls="collapse-paths-' + len + dir + '"><h3 class="display-7">' +
    '<a href="#" class="stmt_toggle"> ' + len + ' Edge Paths (<div id="subject-placeholder-' +
    len + dir + '" class="placeholder subject-placeholder">A</div>&rarr;' + intermArrows +
    '<div id="object-placeholder-' + len + dir + '" class="placeholder object-placeholder">B</div>)</a><span ' +
    'id="npaths-' + len + dir + '" class="badge badge-primary badge-pill float-right path-count">Paths: ' +
    '0</span></h3></div><div id="collapse-paths-' + len + dir + '" class="collapse show"><div class="card-body">' +
    '<table class="table"><thead class="table-head"><th>Path</th><th>Weight</th><th>Support</th></thead><tbody ' +
    'class="table-body" id="query-results-' + len + dir + '"></tbody></table></div></div>';
  return newCard;
}

function generateCommonParents(linkURL) {
  let newCard = document.createElement('div');
  newCard.className = 'card';
  newCard.innerHTML = '<div class="card-header" data-toggle="collapse" ' +
    'data-target="#collapse-paths-cp" aria-expanded="true" aria-controls="collapse-paths-cp">' +
    '<h3 class="display-7"><a href="#" class="stmt_toggle">Complexes and Families (<span id="subject-placeholder-cp" ' +
    'class="placeholder subject-placeholder">A</span> - <span id="object-placeholder-cp" class="placeholder ' +
    'object-placeholder">B</span>)</a><span id="npaths-cp" class="badge badge-primary badge-pill float-right ' +
    'path-count">Entities: 0</span></h3></div><div id="collapse-paths-cp" class="collapse show">' +
    '<div class="card-body"><table class="table"><thead' +
    ` class="table-head"><th id="linkout0"><a title="Right click and copy link to link to this result" class="stmt_toggle" href="${linkURL}#linkout0">Family/Complex</a></th>` +
    '<th>Namespace</th></thead><tbody class="table-body" id="query-results-cp"></tbody></table></div></div>';

  return newCard;
}

function generateCommonTargets() {
  let newCard = document.createElement('div');
  newCard.className = 'card';
  newCard.innerHTML = '<div class="card-header" data-toggle="collapse" data-target="#collapse-common-targets" ' +
    'aria-expanded="true" aria-controls="collapse-common-targets"><h3 class="display-7">' +
    '<a href="#" class="stmt_toggle">Common Targets (<span id="subject-placeholder-ct" ' +
    'class="placeholder subject-placeholder">A</span>&rarr;Z&larr;<span id="object-placeholder-ct" '+
    'class="placeholder object-placeholder">B</span>)</a><span id="common-targets" class="badge badge-primary ' +
    'badge-pill float-right path-count">Targets: 0</span></h3></div><div id="collapse-common-targets" ' +
    'class="collapse show"><div class="card-body"><table class="table"><thead class="table-head"><th>Target ' +
    '(Z)</th><th>Support</th></thead><tbody class="table-body" id="query-results-common-targets"></tbody>' +
    '</table></div></div>';

  return newCard;
}

function generateSharedRegulators() {
  let newCard = document.createElement('div');
  newCard.className = 'card';
  newCard.innerHTML = '<div class="card-header" data-toggle="collapse" data-target="#collapse-shared-regulators" ' +
    'aria-expanded="true" aria-controls="collapse-shared-regulators"><h3 class="display-7">' +
    '<a href="#" class="stmt_toggle"> Shared Regulators (<span id="subject-placeholder-sr" ' +
    'class="placeholder subject-placeholder">A</span>&larr;Z&rarr;<span id="object-placeholder-sr" '+
    'class="placeholder object-placeholder">B</span>)</a><span id="shared-regulators-span" class="badge badge-primary ' +
    'badge-pill float-right path-count">Regulators: 0</span></h3></div><div id="collapse-shared-regulators" ' +
    'class="collapse show"><div class="card-body"><table class="table"><thead class="table-head"><th>Target ' +
    '(Z)</th><th>Support</th></thead><tbody class="table-body" id="query-results-shared-regulators"></tbody>' +
    '</table></div></div>';

  return newCard;
}
