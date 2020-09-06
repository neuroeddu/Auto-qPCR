// Add the following code if you want the name of the file appear on select
document.querySelector(".custom-file-input").addEventListener("change",function(e){
	$(this).next(".custom-file-label").html(e.target.files.length + " file(s) uploaded");
});

// Display and hide file information fields
document.getElementById("info_form").style.display = "none";

document.getElementById("y_fileinfo").onclick = function() {
	document.getElementById("info_form").style.display = "block";
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
	document.getElementById("quantity").value = "";
	document.getElementById("opt_gcol").checked = true;
	document.getElementById("opt_glist").checked = false;
	document.getElementById("glist").value = "";
	document.getElementById("gcol").value = "";
	document.getElementById("y_rm").checked = false;
	document.getElementById("n_rm").checked = true;
	document.getElementById("y_nd").checked = true;
	document.getElementById("n_nd").checked = false;
	document.getElementById("stats_form").style.display = "none";
};

// display or hide one-way or two-way ANOVA
document.getElementById("form-twoway").style.display = "none";
document.getElementById("opt_gcol").setAttribute("checked","");

document.getElementById("twoway").onclick = function() {
	document.getElementById("opt_gcol").checked = false;
	document.getElementById("opt_glist").checked = false;
	document.getElementById("tw_glist").checked = false;
	document.getElementById("tw_gcol").checked = true;
	document.getElementById("gcol1").removeAttribute("disabled");
	document.getElementById("gcol2").removeAttribute("disabled");
	document.getElementById("colname1").setAttribute("disabled","");
	document.getElementById("colname1").value = "";
	document.getElementById("colname2").setAttribute("disabled","");
	document.getElementById("colname2").value = "";
	document.getElementById("glist1").setAttribute("disabled","");
	document.getElementById("glist1").value = "";
	document.getElementById("glist2").setAttribute("disabled","");
	document.getElementById("glist2").value = "";
	document.getElementById("form-oneway").style.display = "none";
	document.getElementById("form-twoway").style.display = "block";
};

document.getElementById("oneway").onclick = function() {
	document.getElementById("tw_gcol").checked = false;
	document.getElementById("opt_glist").checked = false;
	document.getElementById("tw_glist").checked = false;
	document.getElementById("opt_gcol").checked = true;
	document.getElementById("gcol").disabled = false;
	document.getElementById("glist").disabled = true;
	document.getElementById("form-oneway").style.display = "block";
	document.getElementById("form-twoway").style.display = "none";
};

// One-way ANOVA: disable group column or group list when the other option is selected
document.getElementById("glist").setAttribute("disabled","");

document.getElementById("opt_gcol").onclick = function() {
	document.getElementById("gcol").disabled = false;
	document.getElementById("glist").disabled = true;
	document.getElementById("glist").value = "";
};

document.getElementById("opt_glist").onclick = function() {
	document.getElementById("gcol").disabled = true;
	document.getElementById("glist").disabled = false;
	document.getElementById("glist").removeAttribute("disabled");
};

// Two-way ANOVA: disable group column or group list when the other option is selected
document.getElementById("colname1").setAttribute("disabled","");
document.getElementById("colname2").setAttribute("disabled","");
document.getElementById("glist1").setAttribute("disabled","");
document.getElementById("glist2").setAttribute("disabled","");

document.getElementById("tw_gcol").onclick = function() {
	document.getElementById("gcol1").disabled = false;
	document.getElementById("gcol2").disabled = false;
	document.getElementById("colname1").disabled = true;
	document.getElementById("colname1").value = "";
	document.getElementById("colname2").disabled = true;
	document.getElementById("colname2").value = "";
	document.getElementById("glist1").disabled = true;
	document.getElementById("glist1").value = "";
	document.getElementById("glist2").disabled = true;
	document.getElementById("glist2").value = "";
};

document.getElementById("tw_glist").onclick = function() {
	document.getElementById("gcol1").disabled = true;
	document.getElementById("gcol1").value = "";
	document.getElementById("gcol2").disabled = true;
	document.getElementById("gcol2").value = "";
	document.getElementById("colname1").disabled = false;
	document.getElementById("colname2").disabled = false;
	document.getElementById("glist1").disabled = false;
	document.getElementById("glist2").disabled = false;
};

// clear form
document.getElementById("clear").onclick = function() {
	location.reload();
};
