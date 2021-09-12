// Add the following code if you want the name of the file appear on select
document.querySelector(".custom-file-input").addEventListener("change",function(e){
	$(this).next(".custom-file-label").html(e.target.files.length + " file(s) uploaded");
});

// Display and hide file information fields
document.getElementById("info_form").style.display = "none";

document.getElementById("y_fileinfo").onclick = function() {
	document.getElementById("info_form").style.display = "block";
};

document.getElementById("n_fileinfo").onclick = function() {
	document.getElementById("genes").value = "";
	document.getElementById("quencher").value = "";
	document.getElementById("task").value = "UNKNOWN";
	document.getElementById("info_form").style.display = "none";
};

// CSample is only enabled and required when relative ddCT or instability is selected
document.getElementById("csample").setAttribute("disabled", "");

document.getElementById("specialStar").style.display = "none";

document.getElementById("absolute").onclick = function() {
	document.getElementById("csample").setAttribute("disabled", "");
	document.getElementById("specialStar").style.display = "none";
};

document.getElementById("relative_dCT").onclick = function() {
	document.getElementById("csample").setAttribute("disabled", "");
	document.getElementById("specialStar").style.display = "none";
};

document.getElementById("relative_ddCT").onclick = function() {
	document.getElementById("csample").setAttribute("required","");
	document.getElementById("csample").removeAttribute("disabled");
	document.getElementById("specialStar").style.display = "inline-block";
};

document.getElementById("instability").onclick = function() {
	document.getElementById("csample").setAttribute("required","");
	document.getElementById("csample").removeAttribute("disabled");
	document.getElementById("specialStar").style.display = "inline-block";
};

// display or hide stats form
document.getElementById("stats_form").style.display = "none";

document.getElementById("y_stats").onclick = function() {
	document.getElementById("stats_form").style.display = "block";
};

document.getElementById("n_stats").onclick = function() {
	document.getElementById("colnames").value = "";
	document.getElementById("quantity").value = "";
	document.getElementById("opt-gcol").checked = true;
	document.getElementById("opt-glist").checked = false;
	document.getElementById("glist").value = "";
	document.getElementById("glist1").value = "";
	document.getElementById("glist2").value = "";
	document.getElementById("gcol").value = "";
	document.getElementById("gcol1").value = "";
	document.getElementById("gcol2").value = "";
	document.getElementById("colname1").value = "";
	document.getElementById("colname2").value = "";
	document.getElementById("y_rm").checked = false;
	document.getElementById("n_rm").checked = true;
	document.getElementById("y_nd").checked = true;
	document.getElementById("n_nd").checked = false;
	document.getElementById("stats_form").style.display = "none";
};

// disable group column or group list when the other option is selected
// display or hide one-way or two-way ANOVA
document.getElementById("opt-gcol").setAttribute("checked","");
document.getElementById("oneway").setAttribute("checked", "");
document.getElementById("col-form").style.display = "inline";
document.getElementById("list-form").style.display = "none";
document.getElementById("group-col2").style.display = "none";
document.getElementById("group-list2").style.display = "none";

document.getElementById("opt-gcol").onclick = function() {
	document.getElementById("col-form").style.display = "block";
	document.getElementById("list-form").style.display = "none";
	document.getElementById("glist").value = "";
	document.getElementById("colname1").value = "";
	document.getElementById("colname2").value = "";
	document.getElementById("glist1").value = "";
	document.getElementById("glist2").value = "";
};

document.getElementById("opt-glist").onclick = function() {
	document.getElementById("col-form").style.display = "none";
	document.getElementById("list-form").style.display = "block";
	document.getElementById("gcol").value = "";
	document.getElementById("gcol1").value = "";
	document.getElementById("gcol2").value = "";
};

document.getElementById("oneway").onclick = function() {
	document.getElementById("group-col2").style.display = "none";
	document.getElementById("group-list2").style.display = "none";
	document.getElementById("group-col").style.display = "block";
	document.getElementById("group-list").style.display = "block";
	document.getElementById("colname1").value = "";
	document.getElementById("colname2").value = "";
	document.getElementById("glist1").value = "";
	document.getElementById("glist2").value = "";
}

document.getElementById("twoway").onclick = function() {
	document.getElementById("group-col").style.display = "none";
	document.getElementById("group-list").style.display = "none";
	document.getElementById("group-col2").style.display = "block";
	document.getElementById("group-list2").style.display = "block";
	document.getElementById("gcol").value = "";
	document.getElementById("gcol1").value = "";
	document.getElementById("gcol2").value = "";
}

//clear form
document.getElementById("clear").onclick = function() {
	location.reload();
};
