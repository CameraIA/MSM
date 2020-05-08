run("Close All");

path = "/Users/dani/Dropbox/prog/Apps_CIP/JSR_2017/msm/examples/data/99_5p5_region/orig/";

run("Image Sequence...", "open="+path+"res/kmeans/tif sort");
rename("kmeans");
run("Find Edges", "stack");

run("Image Sequence...", "open="+path+"/res/srm/tif sort");
rename("srm");
run("Find Edges", "stack");

/*
run("Image Sequence...", "open="+path+"/res/pmrf/tif sort");
rename("pmrf");
run("Find Edges", "stack");
*/

//Create one stack
//run("Merge Channels...", "c1=kmeans c2=srm c3=pmrf create keep ignore");
run("Merge Channels...", "c1=kmeans c2=srm  create keep ignore");
run("Flatten");

print("Red=kmeans, Green=srm");
setColor(0, 0, 255);
setFont("SansSerif", 20, "bold");
drawString("Red=kmeans, Green=srm",20,20);
//setFont("SansSerif", 40, "bold"); drawString("CAMERA",20,100);
run("Animated Gif ... ", "name=Composite set_global_lookup_table_options=[Do not use] optional=[] image=[No Disposal] set=400 number=0 transparency=[No Transparency] red=0 green=0 blue=0 index=0 filename=/Users/dani/Dropbox/aqui/Math/CAMERA/Dula/GE_fibers/2017_09_anjali_complain/target/analysis/Overlay.gif");
