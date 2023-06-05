// Globals
var gList = [];
var sList = [];
var files = [];

document.addEventListener("DOMContentLoaded", () => {
  // Select Centrality Measure
  const measure = document.getElementById("measure");
  const local = document.getElementById("local");
  const global = document.getElementById("global");
  const query = document.getElementById("query");
  const help = document.getElementById("shelp");
  const globalAdd = document.getElementById("global-add");
  const tableDiv = document.getElementById("tableDiv");
  const invTissueDiv = document.getElementById("invalid-tissue");
  const invFileDiv = document.getElementById("invalid-file");

  // Global table data
  var tBody = document.getElementById("tBody");
  var gTissue = document.getElementById("global-tissue");
  var gDataset = document.getElementById("global-dataset");
  var gFile = document.getElementById("global-file");
  var query_tissue = document.getElementById('query-tissue');
  var query_file = document.getElementById('query-file');

  measure.addEventListener("change", function () {
    let help_content = "";
    switch (measure.value) {
      case "global":
        local.classList.add("d-none");
        global.classList.remove("d-none");
        query.classList.add("d-none");
        help_content = `<h5><i class="fa fa-circle-info"></i> How to use</h5>
        <ul class="list-group list-group-flush bg-transparent">
            <li class="list-group-item">
              <strong>Step 1:</strong> Add tissues to the <strong>Tissue List</strong> to create the
              network by following the following steps:
              <ol class="list-group list-group-numbered">
                  <li class="list-group-item">Select or type <strong>Tissue Name.</strong></li>
                  <li class="list-group-item">Select dataset (default:<strong>GTex</strong>) from <strong>Select Dataset</strong> dropdown. If you want to use your own gene expression data, then select <strong>Others</strong> option.</li>
                  <li class="list-group-item">If you selected <strong>Others</strong> option from <strong>Dataset</strong> dropdown, Then select a compressed (<strong>.zip</strong>, <strong>.gz</strong>) csv file.</li>
              </ol>
            </li>
            <li class="list-group-item">
                <strong>Step 2:</strong> Select the Centrality Measure from <strong> Select Centrality
                    Measure</strong> dropdown.
            </li>
            <li class="list-group-item">
                <strong>Step 3:</strong> Click the <strong>Get Ranking </strong> button.
            </li>
            <li class="list-group-item">
                <strong>Step 4:</strong> Please wait for <strong>3-5 min(s)</strong>.
            </li>
            <li class="list-group-item">
                <strong>Step 5:</strong> Scroll down to see the results.
            </li>
            <li class="list-group-item">
                <strong>Step 6:</strong> You can filter the results using <strong>Search</strong> option.
            </li>
        </ul>
        <div class="container-fluid mt-5">
            <h6><i class="fa fa-download"></i> Sample files used for <strong>Example Run</strong></h6>
            <ul class="list-group list-group-flush bg-transparent">
                <li class="list-group-item">
                    <a href="static/downloads/Pancreas.zip" download><i class="fa fa-download"></i>
                        Pancreas.zip</a> (Compressed csv file containing gene expression data for Pancreas
                    Tissue)
                </li>
                <li class="list-group-item">
                    <a href="static/downloads/Muscle.zip" download><i class="fa fa-download"></i>
                        Muscle.zip</a> (Compressed csv file containing gene expression data for Muscle
                    Tissue)
                </li>
            </ul>
        </div>`;
        break;
      case "query":
        local.classList.add("d-none");
        global.classList.add("d-none");
        query.classList.remove("d-none");
        help_content = `<h5><i class="fa fa-circle-info"></i> How to use</h5>
          <ul class="list-group list-group-flush bg-transparent">
            <li class="list-group-item">
              <strong>Step 1:</strong> Add tissues to the <strong>Tissue List</strong> to create the
              network by following the following steps:
              <ol class="list-group list-group-numbered">
                  <li class="list-group-item">Select or type <strong>Tissue Name.</strong></li>
                  <li class="list-group-item">Select dataset (default:<strong>GTex</strong>) from <strong>Select Dataset</strong> dropdown. If you want to use your own gene expression data, then select <strong>Others</strong> option.</li>
                  <li class="list-group-item">If you selected <strong>Others</strong> option from <strong>Dataset</strong> dropdown, Then select a compressed (<strong>.zip</strong>, <strong>.gz</strong>) csv file.</li>
              </ol>
            </li>
            <li class="list-group-item">
                <strong>Step 2:</strong> Select the Centrality Measure from <strong> Select Centrality
                    Measure</strong> dropdown.
            </li>
            <li class="list-group-item">
              <strong>Step 3:</strong> Select the Query Tissue from <strong>Query Tissue Name</strong> dropdown.
            </li>
            <li class="list-group-item">
              <strong>Step 4:</strong> Select the csv file with gene list (<string>Gene Query Set</string>) for the <strong>Query Tissue</strong>.
            </li>
            <li class="list-group-item">
                <strong>Step 5:</strong> Click the <strong>Get Ranking </strong> button.
            </li>
            <li class="list-group-item">
                <strong>Step 6:</strong> Please wait for <strong>3-5 min(s)</strong>.
            </li>
            <li class="list-group-item">
                <strong>Step 7:</strong> Scroll down to see the results.
            </li>
            <li class="list-group-item">
                <strong>Step 8:</strong> You can filter the results using <strong>Search</strong> option.
            </li>
        </ul>
        <div class="container-fluid mt-5">
            <h6><i class="fa fa-download"></i> Sample files used for <strong>Example Run</strong></h6>
            <ul class="list-group list-group-flush bg-transparent">
                <li class="list-group-item">
                    <a href="static/downloads/Pancreas.zip" download><i class="fa fa-download"></i>
                        Pancreas.zip</a> (Compressed csv file containing gene expression data for Pancreas
                    Tissue)
                </li>
                <li class="list-group-item">
                    <a href="static/downloads/Muscle.zip" download><i class="fa fa-download"></i>
                        Muscle.zip</a> (Compressed csv file containing gene expression data for Muscle
                    Tissue)
                </li>
                <li class="list-group-item">
                    <a href="static/downloads/Pancreas.csv" download><i class="fa fa-download"></i>
                    Pancreas.csv</a> (CSV file containing a list of genes from Pancreas Tissue)
                </li>
            </ul>
        </div>`;
        break;
      default:
        local.classList.remove("d-none");
        global.classList.add("d-none");
        query.classList.add("d-none");
        help_content = `<h5><i class="fa fa-circle-info"></i> How to use</h5>
          <ul class="list-group list-group-flush bg-transparent">
            <li class="list-group-item">
              <strong>Step 1:</strong> Add tissues to the <strong>Tissue List</strong> to create the
              network by following the following steps:
              <ol class="list-group list-group-numbered">
                  <li class="list-group-item">Select or type <strong>Tissue Name.</strong></li>
                  <li class="list-group-item">Select dataset (default:<strong>GTex</strong>) from <strong>Select Dataset</strong> dropdown. If you want to use your own gene expression data, then select <strong>Others</strong> option.</li>
                  <li class="list-group-item">If you selected <strong>Others</strong> option from <strong>Dataset</strong> dropdown, Then select a compressed (<strong>.zip</strong>, <strong>.gz</strong>) csv file.</li>
              </ol>
            </li>
            <li class="list-group-item">
                <strong>Step 2:</strong> Select the Centrality Measure from <strong> Select Centrality
                    Measure</strong> dropdown.
            </li>
            <li class="list-group-item">
                <strong>Step 3:</strong> Click the <strong>Get Ranking </strong> button.
            </li>
            <li class="list-group-item">
                <strong>Step 4:</strong> Please wait for <strong>3-5 min(s)</strong>.
            </li>
            <li class="list-group-item">
                <strong>Step 5:</strong> Scroll down to see the results.
            </li>
            <li class="list-group-item">
                <strong>Step 6:</strong> You can filter the results using <strong>Search</strong> option.
            </li>
        </ul>
        <div class="container-fluid mt-5">
            <h6><i class="fa fa-download"></i> Sample files used for <strong>Example Run</strong></h6>
            <ul class="list-group list-group-flush bg-transparent">
                <li class="list-group-item">
                    <a href="static/downloads/Pancreas.zip" download><i class="fa fa-download"></i>
                        Pancreas.zip</a> (Compressed csv file containing gene expression data for Pancreas
                    Tissue)
                </li>
                <li class="list-group-item">
                    <a href="static/downloads/Muscle.zip" download><i class="fa fa-download"></i>
                        Muscle.zip</a> (Compressed csv file containing gene expression data for Muscle
                    Tissue)
                </li>
            </ul>
        </div>`;
    }
    help.innerHTML = help_content;
  });

  var form_local = document.getElementById("form-local");
  var form_global = document.getElementById("form-global");
  var form_query = document.getElementById("form-query");
  var outputDiv = document.getElementById("result");

  var runSampleButton = document.getElementById("example-run");

  var message = "Please wait! It may take 2 minutes/Tissue  ...";

  runSampleButton.onclick = function (event) {
    loadSampleData();
    var formData = new FormData();
    outputDiv.innerHTML =
      `<img style="display:block;margin-left: auto;margin-right: auto;" height="150px" src="static/img/dna.gif"/><br/>  ${message}`; 
    var formData = new FormData();
    if (measure.value == "global") {
      // var tissue = document.getElementById("global-tissue");
      formData.append("measure", "global");
    } else if (measure.value == "query") {
      query_tissue.selectedIndex = 0;
      formData.append("measure", "query");
      formData.append("tissue-index", 0);
    } else {
      formData.append("measure", "local");
    }
    formData.append('tissues', JSON.stringify(gList));

    var xhr = new XMLHttpRequest();
    xhr.open("POST", "/example_run", true);
    xhr.onload = function () {
      if (xhr.status === 200) {
        outputDiv.innerHTML = xhr.responseText;
        $("#ranking").DataTable({
          columnDefs: [{ className: "dt-center", targets: "_all" }],
          order: [3, "asc"],
        });
      } else {
        outputDiv.innerHTML = "An error occurred during the upload. Try again.";
      }
    };

    xhr.send(formData);
  };

  form_local.onsubmit = function (event) {
    event.preventDefault();

    let pData = prepareData();
    let gList = pData.data;
    let files = pData.files.map(f=>f.file);

    if (gList.length <= 0) {
      $("#error").removeClass("d-none");
      setTimeout(function () {
        $("#error").addClass("d-none");
      }, 5000);
      return;
    }

    outputDiv.innerHTML =
      `<img style="display:block;margin-left: auto;margin-right: auto;" height="150px" src="static/img/dna.gif"/><br/>  ${message}`; 

    var formData = new FormData();
    formData.append("measure", "local");
    formData.append('tissues', JSON.stringify(gList));
    for (let i = 0; i < files.length; i++) {
      if (!(files[i].type.match(".zip") || !files[i].type.match(".gz"))) {
        outputDiv.innerHTML = "Compressed csv file is required!.";
        return;
      }
      if (files[i].size >= 52428800) {
        outputDiv.innerHTML =
          "You cannot upload this file because its size exceeds the maximum limit of 50 MB.";
        return;
      }
      formData.append('files', files[i]);
    }

    var xhr = new XMLHttpRequest();
    xhr.open("POST", "/get_centrality", true);

    // Set up a handler for when the task for the request is complete.
    xhr.onload = function () {
      if (xhr.status === 200) {
        outputDiv.innerHTML = xhr.responseText;
        $("#ranking").DataTable({
          columnDefs: [{ className: "dt-center", targets: "_all" }],
          order: [3, "asc"],
        });
      } else {
        outputDiv.innerHTML = "An error occurred! Please try again.";
      }
    };

    // Send the data.
    xhr.send(formData);
  };

  form_global.onsubmit = function (event) {
    event.preventDefault();

    let pData = prepareData();
    let gList = pData.data;
    let files = pData.files.map(f=>f.file);

    if (gList.length <= 0) {
      $("#error").removeClass("d-none");
      setTimeout(function () {
        $("#error").addClass("d-none");
      }, 5000);
      return;
    }

    outputDiv.innerHTML =
      `<img style="display:block;margin-left: auto;margin-right: auto;" height="150px" src="static/img/dna.gif"/><br/>  ${message}`; 

    var formData = new FormData();
    formData.append("measure", "global");
    formData.append('tissues', JSON.stringify(gList));
    for (let i = 0; i < files.length; i++) {
      if (!(files[i].type.match(".zip") || !files[i].type.match(".gz"))) {
        outputDiv.innerHTML = "Compressed csv file is required!.";
        return;
      }
      if (files[i].size >= 52428800) {
        outputDiv.innerHTML =
          "You cannot upload this file because its size exceeds the maximum limit of 50 MB.";
        return;
      }
      formData.append('files', files[i]);
    }

    var xhr = new XMLHttpRequest();
    xhr.open("POST", "/get_centrality", true);

    // Set up a handler for when the task for the request is complete.
    xhr.onload = function () {
      if (xhr.status === 200) {
        outputDiv.innerHTML = xhr.responseText;
        $("#ranking").DataTable({
          columnDefs: [{ className: "dt-center", targets: "_all" }],
          order: [3, "asc"],
        });
      } else {
        outputDiv.innerHTML = "An error occurred! Please try again.";
      }
    };

    // Send the data.
    xhr.send(formData);
  };

  form_query.onsubmit = function (event) {
    event.preventDefault();

    let pData = prepareData();
    let gList = pData.data;
    let files = pData.files.map(f=>f.file);

    if (gList.length <= 0) {
      $("#error").removeClass("d-none");
      setTimeout(function () {
        $("#error").addClass("d-none");
      }, 5000);
      return;
    }

    outputDiv.innerHTML =
      `<img style="display:block;margin-left: auto;margin-right: auto;" height="150px" src="static/img/dna.gif"/><br/>  ${message}`; 

    if (!(query_file.files[0].type.match(".csv"))) {
      outputDiv.innerHTML = "CSV file is required!.";
      return;
    }
    if (query_file.files[0].size >= 52428800) {
      outputDiv.innerHTML =
        "You cannot upload Gene Query Set file because its size exceeds the maximum limit of 50 MB.";
      return;
    }
    var formData = new FormData();
    formData.append("measure", "query");
    formData.append("query-tissue", query_tissue.value);
    formData.append("tissue-index", query_tissue.selectedIndex);
    formData.append("query-file", query_file.files[0]);
    formData.append('tissues', JSON.stringify(gList));
    for (let i = 0; i < files.length; i++) {
      if (!(files[i].type.match(".zip") || !files[i].type.match(".gz"))) {
        outputDiv.innerHTML = "Compressed csv file is required!.";
        return;
      }
      if (files[i].size >= 52428800) {
        outputDiv.innerHTML =
          "You cannot upload this file because its size exceeds the maximum limit of 50 MB.";
        return;
      }
      formData.append('files', files[i]);
    }

    var xhr = new XMLHttpRequest();
    xhr.open("POST", "/get_centrality", true);

    // Set up a handler for when the task for the request is complete.
    xhr.onload = function () {
      if (xhr.status === 200) {
        outputDiv.innerHTML = xhr.responseText;
        $("#ranking").DataTable({
          columnDefs: [{ className: "dt-center", targets: "_all" }],
          order: [3, "asc"],
        });
      } else {
        outputDiv.innerHTML = "An error occurred! Please try again.";
      }
    };

    // Send the data.
    xhr.send(formData);
  };

  gDataset.onchange = function () {
    if (gDataset.value == "GTex") {
      gFile.required = false;
      gFile.disabled = true;
    } else {
      gFile.required = true;
      gFile.disabled = false;
    }
  };

  gTissue.onchange = function (event) {
    if (gTissue.value.trim() == "") {
      invTissueDiv.style.display = "block";
      gTissue.style.borderColor = "#dc3545";
    } else {
      invTissueDiv.style.display = "none";
      gTissue.style.borderColor = "#ced4da";
    }
  };

  gFile.onchange = function (event) {
    if (gFile.files.length == 0 && gDataset.value != "GTex") {
      invFileDiv.style.display = "block";
      gFile.style.borderColor = "#dc3545";
    } else {
      invFileDiv.style.display = "none";
      gFile.style.borderColor = "#ced4da";
    }
  };

  globalAdd.onclick = function (event) {
    if (gTissue.value.trim() == "") {
      invTissueDiv.style.display = "block";
      gTissue.style.borderColor = "#dc3545";
      return;
    }
    if (gDataset.value != "GTex" && gFile.files.length == 0) {
      invFileDiv.style.display = "block";
      gFile.style.borderColor = "#dc3545";
      return;
    }
    const obj = {
      name: gTissue.value,
      dataset: gDataset.value,
      file: gFile.files.length > 0 ? gFile.files[0].name : "-",
      checked: true
    };
    const isPresent = gList.some(
      (e) =>
        e.name === obj.name && e.dataset === obj.dataset && e.file === obj.file
    );
    if (isPresent) {
      $("#info").removeClass("d-none");
      setTimeout(function () {
        $("#info").addClass("d-none");
      }, 5000);
    } else {
      gList.push(obj);
      if (gFile.files.length > 0) {
        files.push({ id: gList.length - 1, file: gFile.files[0] });
      }
      reloadTable(tBody, gList);
      if (tableDiv.classList.contains("d-none")) {
        tableDiv.classList.remove("d-none");
      }
      gTissue.value = "";
      gDataset.value = "GTex";
      gFile.value = "";
      gFile.disabled = true;
    }
  };
});

function reloadTable(tBody, rows) {
  while (tBody.firstChild) {
    tBody.removeChild(tBody.firstChild);
  }

  var query_tissue = document.getElementById('query-tissue');
  while (query_tissue.firstChild) {
    query_tissue.removeChild(query_tissue.firstChild);
  }

  for (var i = 0; i < rows.length; i++) {
    // Adding Tissue data to table
    var row = document.createElement("tr");
    // Adding checkbox
    var checkboxCell = document.createElement("td");
    checkboxCell.innerHTML = `<input type="checkbox" onchange="checkedChanged(this)" id="gcbx${i}" checked/>`;
    row.appendChild(checkboxCell);
    // Adding Tissue name
    var nameCell = document.createElement("td");
    nameCell.textContent = capitalize(rows[i].name);
    row.appendChild(nameCell);
    // Adding dataset name
    var datasetCell = document.createElement("td");
    datasetCell.textContent = capitalize(rows[i].dataset);
    row.appendChild(datasetCell);
    // Adding file name
    var fileCell = document.createElement("td");
    fileCell.textContent = rows[i].file;
    row.appendChild(fileCell);
    // Adding delete button
    var deleteCell = document.createElement("td");
    deleteCell.innerHTML =
      `<button type="button" onclick = "deleteBtnClick(this)" id="gdb${i}" style="border: none; border-radius: 100%;"><i class="fa fa-trash" aria-hidden="true"></i></button>`;
    row.appendChild(deleteCell);
    // Adding ro to table
    tBody.appendChild(row);

    // Create options for Query Tissue
    var option = document.createElement("option");
    option.value = rows[i].name;
    option.text = capitalize(rows[i].name);
    query_tissue.appendChild(option);
  }
}

function deleteBtnClick(sender) {
  let id = sender.id;
  let index = id[id.length - 1];
  let symbol = id.substring(0, id.length - 1);
  if (symbol == "gdb") {
    removeRow(tBody, gList, index, tableDiv);
  } else {
    // removeRow(tBody, gList, index, tableDiv);
  }
}

function checkedChanged(sender) {
  let id = sender.id;
  let index = id[id.length - 1];
  let symbol = id.substring(0, id.length - 1);
  if (symbol == "gcbx") {
    gList[index].checked = sender.checked;
  }
}

function removeRow(body, rows, index, tableDiv) {
  rows.splice(index, 1);

  reloadTable(body, rows);
  // if (body.children.length == 0) {
  //   tableDiv.classList.add("d-none");
  // }
}

function capitalize(str) {
  return str.charAt(0).toUpperCase() + str.slice(1);
}

function loadSampleData() {
  gList = [];
  files = [];
  gList.push({
    name: 'Pancreas',
    dataset: 'GTex',
    file: "-",
    checked: true
  });
  gList.push({
    name: 'Muscle',
    dataset: 'GTex',
    file: "-",
    checked: true
  });
  reloadTable(tBody, gList);
  if (tableDiv.classList.contains("d-none")) {
    tableDiv.classList.remove("d-none");
  }
}

function prepareData() {
  let sList = [];
  let fList = Array.from(files);
  for (let i = 0; i < gList.length; i++) {
    const e = gList[i];
    if (e.checked) {
      sList.push(e);
    }
    else {
      fList = fList.filter(f => f.id != i);
    }
  }
  return { data: sList, files: fList };
}
