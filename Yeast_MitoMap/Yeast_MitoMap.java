/*
Yeast_MitoMap - by Richard Butler, Gurdon Institute Imaging Facility, University of Cambridge
Copyright 2013, 2014 Richard Butler

Yeast_MitoMap is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Yeast_MitoMap is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Yeast_MitoMap.  If not, see <http://www.gnu.org/licenses/>.
 */

import ij.*;
import ij.plugin.*;
import ij.gui.*;
import ij.measure.*;

import ij3d.*;

import java.awt.Font;
import java.util.ArrayList;
import java.util.Arrays;

import javax.swing.*;
import javax.vecmath.Color3f;

import java.awt.*;
import java.awt.event.*;


public class Yeast_MitoMap implements PlugIn,ActionListener{
	private ImagePlus imp, map;
	private JFrame gui;
	private JCheckBox check3D;
	private boolean do3D = false;
	private Roi userRoi;
	private Overlay ol;
	private int W,H,C,Z,T, roiX,roiY,roiW,roiH, objC;
	private Calibration cal;
	private String title,unit;
	private double[] MVarea = new double[16];
	private static final Font labelFont = new Font("Helvetica",Font.BOLD,12);
	private final Image icon = Toolkit.getDefaultToolkit().getImage(getClass().getResource("logo_icon.gif"));

	public void run(String arg){
		if(WindowManager.getImageCount()==0){IJ.error("No Image","No Images are open.");return;}
		imp = WindowManager.getCurrentImage();
		title = imp.getTitle();
		cal = imp.getCalibration();
		W = imp.getWidth();	H = imp.getHeight();
		C = imp.getNChannels();	Z = imp.getNSlices();	T = imp.getNFrames();
		unit = imp.getCalibration().getUnit();
		objC = C;	//analyse last channel in stack
		if(unit.matches("[Mm]icrons?")){unit="\u00B5m";}

		//area contributions of voxel classes
		//factors from Mullikin and Verbeek (2001), see G.Windreich et al. (2003) for open access review
		//extra classes added to deal with different areas for horizontal and vertical faces
		MVarea[1] = cal.pixelWidth*cal.pixelDepth*0.894d;
		MVarea[2] = cal.pixelWidth*cal.pixelDepth*1.3409d;
		MVarea[3] = (((2d*(cal.pixelWidth*cal.pixelDepth))+(cal.pixelWidth*cal.pixelHeight))/3d)*1.5879d;
		MVarea[4] = cal.pixelWidth*cal.pixelDepth*2d;
		MVarea[5] = (((2d*cal.pixelWidth*cal.pixelDepth)+(2d*cal.pixelWidth*cal.pixelHeight))/4d)*(8d/3d);
		MVarea[6] = (((3d*(cal.pixelWidth*cal.pixelDepth))+(2d*cal.pixelWidth*cal.pixelHeight))/3d)*(10d/3d);
		MVarea[7] = 2d*cal.pixelWidth*cal.pixelDepth;
		MVarea[8] = 4d*cal.pixelWidth*cal.pixelDepth;
		MVarea[9] = (4d*cal.pixelWidth*cal.pixelDepth)+(2d*cal.pixelWidth*cal.pixelHeight);
		MVarea[10] = cal.pixelWidth*cal.pixelDepth;
		MVarea[11] = (((cal.pixelWidth*cal.pixelDepth)+((cal.pixelWidth*cal.pixelHeight)))/2d)*1.3409d;
		MVarea[12] = (((3d*cal.pixelWidth*cal.pixelDepth)+(cal.pixelWidth*cal.pixelHeight))/4d)*(8d/3d);
		MVarea[13] = (((4d*(cal.pixelWidth*cal.pixelDepth))+(cal.pixelWidth*cal.pixelHeight))/3d)*(10d/3d);
		MVarea[14] = 2d*cal.pixelWidth*cal.pixelHeight;
		MVarea[15] = ((((cal.pixelWidth*cal.pixelDepth))+(2d*cal.pixelWidth*cal.pixelHeight))/3d)*1.5879d;

		gui = new JFrame("Yeast_MitoMap : "+title);
		gui.setIconImage(icon);
		gui.setLayout(new BorderLayout());;
		JPanel checkPanel = new JPanel();
		check3D = new JCheckBox("Show in 3D Viewer", do3D);
		checkPanel.add(check3D);
		gui.add(checkPanel,BorderLayout.CENTER);
		JPanel buttons = new JPanel();
		JButton goButton = new JButton("Analyse");
		goButton.addActionListener(this);
		JButton closeButton = new JButton("Close");
		closeButton.addActionListener(this);
		buttons.add(goButton);
		buttons.add(closeButton);
		gui.add(buttons,BorderLayout.SOUTH);
		gui.pack();
		gui.setLocation((int)Math.round(Toolkit.getDefaultToolkit().getScreenSize().getWidth()/2),200);
		gui.setVisible(true);
	}

	private void mapMitochondria(){
		try{
			//create object mask
			Prefs.blackBackground = true;
			map = new Duplicator().run(imp,objC,objC,1,Z,1,1);
			map.show();
			IJ.run(map, "16-bit", "");
			//IJ.run(map, "Subtract Background...", "rolling=10 stack");
			IJ.setAutoThreshold(map, "Otsu dark stack");
			IJ.run(map, "Convert to Mask", "method=Otsu background=Dark black");
			//IJ.run(map, "Remove Outliers...", "radius=1 threshold=0 which=Bright stack");
			map.setRoi(userRoi);
			map.getRoi().setLocation(0,0);
			IJ.run(map, "Clear Outside", "stack");
			IJ.run(map, "Select None", "");
			//if(1==1){map.show();return;}	//debugging
			if(imp.getOverlay()==null){
				ol = new Overlay();
			}
			else{
				ol = imp.getOverlay();
			}
			userRoi.setPosition(0,0,0);
			userRoi.setStrokeColor(Color.MAGENTA);
			ol.add(userRoi);

			FindConnectedRegions fcr = new FindConnectedRegions();
			fcr.run(map,true,true,false,false,false,false,false,0d,50,-1,false);

			int[] id = WindowManager.getIDList();
			ResultsTable rt = new ResultsTable();
			rt.showRowNumbers(false);
			int ref = -1;
			for(int i=0;i<id.length;i++){
				ImagePlus obj = WindowManager.getImage(id[i]);
				String objTitle = obj.getTitle();
				int objW = obj.getWidth();
				int objH = obj.getHeight();
				if(!objTitle.matches("^Region of value .*")){continue;}
				ref++;
				Prefs.blackBackground = true;
				ArrayList<Point3b> vox = new ArrayList<Point3b>();
				ArrayList<Point3b> voxCal = new ArrayList<Point3b>();
				double centroidX = 0d; double centroidY = 0d; double centroidZ = 0d;
				double sum = 0d; int n = 0;
				for(int z=1;z<=Z;z++){
					obj.setSlice(z);
					imp.setPositionWithoutUpdate(objC,z,1);
					Prefs.blackBackground = true;
					for(int y=0;y<objH;y++){
						for(int x=0;x<objW;x++){
							if(obj.getPixel(x,y)[0]==255){
								vox.add(new Point3b(x,y,z));
								voxCal.add(new Point3b(x*cal.pixelWidth,y*cal.pixelHeight,z*cal.pixelDepth));
								centroidX += x;
								centroidY += y;
								centroidZ += z;
								sum += Float.intBitsToFloat(imp.getPixel(roiX+x,roiY+y)[0]);
								n++;
							}
						}
					}
				}
				//IJ.log(mean+" / "+n);
				double mean = sum/n;
				//IJ.log("= "+mean);
				centroidX = centroidX/vox.size();
				centroidY = centroidY/vox.size();
				centroidZ = centroidZ/vox.size();
				Point3b centroid = new Point3b(centroidX,centroidY,centroidZ);
				Point3b centroidCal = new Point3b(centroidX*cal.pixelWidth,centroidY*cal.pixelHeight,centroidZ*cal.pixelDepth);

				ArrayList<Point3b> surface = new ArrayList<Point3b>();
				ArrayList<Point3b> surfaceCal = new ArrayList<Point3b>();
				double surfaceArea = 0d;
				for(int a=0;a<vox.size();a++){
					int Hcount = 2;	int Vcount = 4;
					for(int b=0;b<vox.size();b++){
						if(a==b){continue;}
						if(vox.get(a).distance(vox.get(b))==1){	//1 is uncalibrated distance to face-connected neighbour
							boolean deep = vox.get(a).z != vox.get(b).z;
							//count horizontal and vertical connected faces for voxel classification
							if(deep){Hcount--;}
							else{Vcount--;}
						}
					}
					if(Hcount+Vcount>0){	//if this voxel has exposed faces
						surface.add(vox.get(a));
						surfaceCal.add(new Point3b(vox.get(a).x*cal.pixelWidth,vox.get(a).y*cal.pixelHeight,vox.get(a).z*cal.pixelDepth));
						int MVclass = -1;
						if	(Hcount==0&&Vcount==1){MVclass=1;}
						else if (Hcount==0&&Vcount==2){MVclass=2;}
						else if (Hcount==1&&Vcount==2){MVclass=3;}
						else if (Hcount==0&&Vcount==3){MVclass=4;}
						else if (Hcount==2&&Vcount==2){MVclass=5;}
						else if (Hcount==2&&Vcount==3){MVclass=6;}
						//else if (Hcount==0&&Vcount==2){MVclass=7;}	//TODO? handle two opposite exposed faces
						else if (Hcount==0&&Vcount==4){MVclass=8;}
						else if (Hcount==2&&Vcount==4){MVclass=9;}
						else if (Hcount==1&&Vcount==0){MVclass=10;}
						else if (Hcount==1&&Vcount==1){MVclass=11;}
						else if (Hcount==1&&Vcount==3){MVclass=12;}
						else if (Hcount==1&&Vcount==4){MVclass=13;}	
						else if (Hcount==2&&Vcount==0){MVclass=14;}
						else if (Hcount==2&&Vcount==1){MVclass=15;}
						else{throw new RuntimeException("Unclassified Voxel : "+Hcount+":"+Vcount);}
						surfaceArea += MVarea[MVclass];
						Roi mark = new Roi(roiX+vox.get(a).x,roiY+vox.get(a).y,1,1);
						mark.setPosition(0,(int)Math.round(vox.get(a).z),0);
						mark.setFillColor(Color.MAGENTA);
						ol.add(mark);
					}
				}

				double meanDtoCentroid = 0d;
				double meanDtoXaxis = 0d;
				double meanDtoYaxis = 0d;
				double meanDtoZaxis = 0d;
				for(int v=0;v<vox.size();v++){
					meanDtoCentroid += voxCal.get(v).distance(centroidCal);
					meanDtoXaxis += voxCal.get(v).distance(new Point3b(voxCal.get(v).x,centroidCal.y,centroidCal.z));
					meanDtoYaxis += voxCal.get(v).distance(new Point3b(centroidCal.x,voxCal.get(v).y,centroidCal.z));
					meanDtoZaxis += voxCal.get(v).distance(new Point3b(centroidCal.x,centroidCal.y,voxCal.get(v).z));
				}
				meanDtoCentroid = meanDtoCentroid/vox.size();
				meanDtoXaxis = meanDtoXaxis/vox.size();
				double varDtoCentroid = 0d;	//variance of distance to the centroid = second central moment
				double varDtoXaxis = 0d;	//second moments taking the centroid as the origin
				double varDtoYaxis = 0d;
				double varDtoZaxis = 0d;
				for(int v=0;v<vox.size();v++){
					varDtoCentroid += Math.pow(voxCal.get(v).distance(centroidCal)-meanDtoCentroid,2);
					varDtoXaxis += Math.pow(voxCal.get(v).distance(new Point3b(voxCal.get(v).x,centroidCal.y,centroidCal.z))-meanDtoXaxis,2);
					varDtoYaxis += Math.pow(voxCal.get(v).distance(new Point3b(centroidCal.x,voxCal.get(v).y,centroidCal.z))-meanDtoYaxis,2);
					varDtoZaxis += Math.pow(voxCal.get(v).distance(new Point3b(centroidCal.x,centroidCal.y,voxCal.get(v).z))-meanDtoZaxis,2);
				}
				varDtoCentroid = varDtoCentroid/vox.size();
				varDtoXaxis = varDtoXaxis/vox.size();
				varDtoYaxis = varDtoYaxis/vox.size();
				varDtoZaxis = varDtoZaxis/vox.size();

				double meanRadius = 0d;
				for(int s=0;s<surfaceCal.size();s++){
					meanRadius += surfaceCal.get(s).distance(centroidCal);
				}
				meanRadius = meanRadius/surfaceCal.size();
				double varRadius = 0d;
				for(int s=0;s<surfaceCal.size();s++){
					varRadius += Math.pow(surfaceCal.get(s).distance(centroidCal)-meanRadius,2);
				}
				varRadius = varRadius/surfaceCal.size();

				double volume = Double.valueOf(vox.size())*(cal.pixelWidth*cal.pixelHeight*cal.pixelDepth);
				double compactnessCP = varDtoCentroid/volume;	//cell profiler method
				double ratio1 = Math.min(varDtoXaxis,varDtoYaxis)/Math.max(varDtoXaxis,varDtoYaxis);
				double ratio2 = Math.min(varDtoXaxis,varDtoZaxis)/Math.max(varDtoXaxis,varDtoZaxis);
				double ratio3 = Math.min(varDtoYaxis,varDtoZaxis)/Math.max(varDtoYaxis,varDtoZaxis);
				double iso =  ratio1+ratio2+ratio3;
				double sphereSA = (Math.pow(Math.PI,1d/3d)*Math.pow(6d*volume,2d/3d));	//surface area of a sphere with object volume	
				double sphereR = Math.sqrt(surfaceArea/(4d*Math.PI));	//radius of a sphere with object surface area
				double sphereV = (4d/3d)*Math.PI*Math.pow(sphereR,3);	//volume of a sphere with object surface area
				/* double sphereR = Math.cbrt((3d*volume)/(4d*Math.PI));//radius of a sphere with object volume
		double sphereV = (4d/3d)*Math.PI*Math.pow(sphereR,3);	//volume of a sphere with object volume
		double sphereSA = 4d*Math.PI*Math.pow(sphereR,2);	//surface area of a sphere with object volume */
				double ipq = volume/sphereV;	//isoperimetric quotient
				double sphericity = sphereSA/surfaceArea;
				double SAVratio = surfaceArea/volume;

				int row = rt.getCounter();
				rt.setValue("Object",row,ref);
				rt.setLabel(title,row);
				rt.setValue("Object X",row,(roiX+centroidX)*cal.pixelWidth);
				rt.setValue("Object Y",row,(roiY+centroidY)*cal.pixelHeight);
				rt.setValue("Object Z",row,centroidZ*cal.pixelDepth);
				rt.setValue("Mean",row,mean);
				rt.setValue("Sum",row,sum);
				rt.setValue("Volume ("+unit+"\u00b3)",row,volume);
				rt.setValue("Surface Area ("+unit+"\u00b2)",row,surfaceArea);
				rt.setValue("Compactness (CP)",row,compactnessCP);
				rt.setValue("Distribution Isotropy",row,iso);
				rt.setValue("Isoperimetric Quotient",row,ipq);
				rt.setValue("Sphericity",row,sphericity);
				rt.setValue("SA:V",row,SAVratio);
				rt.setValue("Radius Variance",row,varRadius);
				String rtTitle = "Yeast_MitoMap : "+title+" ("+IJ.d2s((roiX+(roiW/2))*cal.pixelWidth,2)+","+IJ.d2s((roiY+(roiH/2))*cal.pixelWidth,2)+")";
				rt.show(rtTitle);
				WindowManager.getFrame(rtTitle).setBounds(new Rectangle(100,100,1600,400));

				for(int z=1;z<=Z;z++){
					obj.setPositionWithoutUpdate(objC,z,1);
					imp.setPositionWithoutUpdate(objC,z,1);
					Prefs.blackBackground = true;
					IJ.run(obj, "Create Selection", "");
					if(obj.getRoi()==null){continue;}
					IJ.run(obj, "Make Inverse", "");
					if(obj.getRoi()==null){continue;}
					Roi roi = obj.getRoi();
					Prefs.blackBackground = true;
					Rectangle bounds = roi.getBounds();
					roi.setLocation(bounds.x+roiX,bounds.y+roiY);
					roi.setPosition(0,z,0);
					roi.setStrokeColor(Color.MAGENTA);
					ol.add(roi);
					IJ.run(imp, "Select None", "");
					IJ.run(map, "Select None", "");
				}

				Roi label = new TextRoi(roiX+centroidX,roiY+centroidY,""+ref,labelFont);
				label.setStrokeColor(Color.CYAN);
				label.setPosition(0,(int)Math.round(centroidZ),0);
				ol.add(label);
				imp.setOverlay(ol);
				obj.changes = false;
				obj.close();
			}
			if(do3D){
				Image3DUniverse univ = new Image3DUniverse();
				Content cont = new Content("");
				boolean[] true3 = {true,true,true};
				cont = univ.addVoltex(map,new Color3f(java.awt.Color.WHITE),title,1,true3,1);
				univ.show();
			}
			map.changes = false;
			map.close();
		}catch(Exception e){IJ.log("map: "+e.toString()+"\n"+Arrays.toString(e.getStackTrace()).replace(",","\n"));}
	}

	public void actionPerformed(ActionEvent event){
		String e = event.getActionCommand();
		if(e=="Analyse"){
			if(imp.getRoi()==null){IJ.error("No Roi","Select an area for analysis.");return;}
			userRoi = imp.getRoi();
			int type = userRoi.getType();
			if((type!=Roi.FREEROI)&&(type!=Roi.OVAL)&&(type!=Roi.POLYGON)&&(type!=Roi.RECTANGLE)&&(type!=Roi.TRACED_ROI)){
				IJ.error("Area Required","Select an area for analysis.");
				return;
			}
			do3D = check3D.isSelected();
			gui.dispose();
			if(type==Roi.FREEROI){IJ.run("Fit Spline", "");userRoi = imp.getRoi();}
			Rectangle userBounds = userRoi.getBounds();
			roiX = userBounds.x;
			roiY = userBounds.y;
			roiW = userBounds.width;
			roiH = userBounds.height;
			mapMitochondria();
		}
		else if(e=="Close"){
			gui.dispose();
		}
	}

}
