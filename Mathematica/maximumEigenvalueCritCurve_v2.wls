#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Input:: *)
(*Solve[{r==b-d,AI/PI==d/b},{b,d}]//FullSimplify*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*beta[x_, \[Eta]_] := 1-1/(\[Eta] + E^(Abs[x-2]^0.5))*)
(*(*beta[x_, \[Eta]_] := \[Eta]*)*)


(* ::Input:: *)
(*beta[2,99.001]*)


(* ::Input:: *)
(*Plot[{beta[x,0.01],beta[x,0.5], beta2[x, 0.4],beta2[x, 0.2]},{x,1,8}, PlotStyle->{Blue,Darker[Blue],Red,Darker[Red]},Method->"DefaultPlotStyle"->Thick,Frame->True,FrameStyle->Directive[Thickness[0.0035],Black],FrameLabel->{"Ploidy", "Mis-segregration Rate (\[Beta])"},LabelStyle->{FontSize->18,Black,FontFamily->"Arial"},PlotRange->{{1,8},Automatic}, AspectRatio->1]*)


(* ::Input:: *)
(*myfun[\[Lambda]_,\[Delta]_?NumberQ,n_?NumberQ, \[Eta]_ ?NumberQ]:=Eigenvalues[Table[If[i==j, \[Lambda](2(1-beta[j, \[Eta]])^i-1)-\[Delta],If[2j>=i>j,\[Lambda] Binomial[j,i-j]beta[j, \[Eta]]^(i-j) (1-beta[j, \[Eta]])^(2j-i),If[j>i,\[Lambda] Binomial[j,j-i]beta[j, \[Eta]]^(j-i) (1-beta[j, \[Eta]])^i,0]]],{i,n},{j,n}],*)
(*1,Method->{"Arnoldi","Criteria"->"RealPart"}][[1]]*)


(* ::Input:: *)
(*critDelta[\[Delta]_?NumberQ, n_?NumberQ]:=\[Eta]/.FindRoot[myfun[1,\[Delta],n,\[Eta]]==0,{\[Eta],0.9999,0,100}];*)


(* ::Input:: *)
(*nvals={3, 5,8,20};*)
(*ploidies={2};*)
(*leg="n = "<>ToString[#]&/@nvals;*)
(*leg2="p = "<>ToString[#]&/@ploidies;*)
(*colorList=ColorData[37,"ColorList"];*)


(* ::Input:: *)
(*critTable=Table[{\[Delta],beta[2,critDelta[\[Delta],#]]},{\[Delta],0.0001,1,0.1}]&/@nvals;*)
(*(*critTable4=Table[{\[Delta],beta[#,critDelta[\[Delta],8]]},{\[Delta],0.0,1,0.1}]&/@ploidies;*)*)
(*critTableeta=Table[{\[Delta],critDelta[\[Delta],8]},{\[Delta],0.0001,1,0.1}];*)
(*Export[NotebookDirectory[]<>"criticalcurve_EtaVsTurnover"<>".csv",critTableeta[[#]],"CSV"]&/@Range[10]*)


(* ::Input:: *)
(*critTableeta*)


(* ::Input:: *)
(*ListPlot[#&/@critTable,*)
(*Joined->True,*)
(*Frame->True,*)
(*Axes->False,*)
(*PlotRange->{All,All},*)
(*(*Epilog\[Rule]points,*)*)
(*PlotStyle->{{Black},{Black, Dashed}, {Black, DotDashed}, {Black, Dotted}},*)
(*FrameLabel->{"\[Mu]/\[Lambda]","Missegregation rate"},*)
(*FrameStyle->Directive[Thickness[0.0035],Black],*)
(*LabelStyle->{FontSize->18,Black,FontFamily->"Arial"},*)
(*PlotLegends->leg,*)
(*ImageSize->500*)
(*]*)


(* ::Input:: *)
(*ListPlot[#&/@critTable4,*)
(*Joined->True,*)
(*Frame->True,*)
(*Axes->False,*)
(*PlotRange->{All,All},*)
(*PlotStyle->{{Black},{Black, Dashed}, {Black, DotDashed}},*)
(*FrameLabel->{"\[Mu]/\[Lambda]","Missegregation rate"},*)
(*FrameStyle->Directive[Thickness[0.0035],Black],*)
(*LabelStyle->{FontSize->18,Black,FontFamily->"Arial"},*)
(*PlotLegends->{"\[Beta](2)","\[Beta](3)","\[Beta](4)"},*)
(*ImageSize->500*)
(*]*)


(* ::Input:: *)
(*ListPlot[#&/@critTableeta,*)
(*Joined->True,*)
(*Frame->True,*)
(*Axes->False,*)
(*PlotRange->{All,All},*)
(*(*Epilog\[Rule]points,*)*)
(*PlotStyle->{{Red},{Blue}, {Blue, Dashed}, {Black, Dotted}},*)
(*FrameLabel->{"\[Mu]/\[Lambda]","Mis-segregation functional form (\[Eta])"},*)
(*FrameStyle->Directive[Thickness[0.0035],Black],*)
(*LabelStyle->{FontSize->18,Black,FontFamily->"Arial"},*)
(*PlotLegends->nvals,*)
(*ImageSize->500,*)
(*AspectRatio->1*)
(*]*)


(* ::Input:: *)
(*ListPlot[nvals,Joined->True,Frame->True,Axes->False,PlotRange->{All,All},PlotStyle->{{GrayLevel[0]},{GrayLevel[0],Dashing[{Small,Small}]},{GrayLevel[0],Dashing[{0,Small}]},{GrayLevel[0],Dashing[{0,Small,Small,Small}]}},FrameLabel->{"\[Mu]/\[Lambda]","Missegregation rate"},FrameStyle->Directive[Thickness[0.0035`],GrayLevel[0]],LabelStyle->{FontSize->18,GrayLevel[0],FontFamily->"Arial"},PlotLegends->leg,ImageSize->500]*)


(* ::Input:: *)
(*With[{n=5},Table[If[i==j, \[Lambda](1-2\[Mu])-\[Delta],If[Abs[i-j]==1,\[Lambda] \[Mu],0]],{i,n},{j,n}]]//MatrixForm*)


(* ::Input:: *)
(*list={{0.96354,0.01},{0.256,0.01},{0.9874,0.01},{0.9834,0.01},{0.9754,0.01},{0.992,0.01},{0.9496,0.01},{0.9344,0.01},{0.6349,0.01}};*)


(* ::Input:: *)
(*points=Flatten[{PointSize[0.015],colorList[[#]],Point[{list[[#,1]],list[[#,2]]}]}&/@Range[Length[list]],1]*)


(* ::Input:: *)
(*ListLogLinearPlot[critTable,*)
(*Joined->True,*)
(*Frame->True,*)
(*Axes->False,*)
(*PlotRange->All,*)
(*Epilog->points,*)
(*FrameLabel->{"\[Mu]/\[Lambda]","Missegregation rate"},*)
(*FrameStyle->Directive[Thickness[0.0035],Black],*)
(*LabelStyle->{FontSize->18,Black,FontFamily->"Arial"},*)
(*ImageSize->500*)
(*]*)


(* ::Input:: *)
(*makeErrors[med_,low_,high_]:=Around[med,{med-low,high-med}];*)
(*makeErrorsReverse[med_,low_,high_]:=Around[med,{high-med,med-low}];*)


(* ::Input:: *)
(*birthList={{0.192,0.015,1.667},{0.207,0.018,0.625},{0.239,0.019,1.111},{0.303,0.179,0.417},{0.122,0.021,0.556},{0.25,0.159,0.323},{0.139,0.024,0.286},{0.122,0.036,0.667},{0.063,0.043,0.4}};*)
(*growthList={{0.007,0.003,0.033},{0.154,0.035,0.347},{0.003,0.0003,0.038},{0.005,0.0004,0.013},{0.003,0.0003,0.016},{0.002,0.0004,0.012},{0.007,0.002,0.014},{0.008,0.005,0.013},{0.023,0.01,0.5}};*)


(* ::Input:: *)
(*ListLogPlot[makeErrors@@#&/@birthList,*)
(*Ticks->{{},Automatic},*)
(*PlotRange->All,*)
(*ImageSize->500]*)


(* ::Input:: *)
(*makeErrors@@#&/@birthList*)


(* ::Input:: *)
(*labels={"Head and neck","Esophogaeal","Colorectal","Rectal","Breast","Cervix","Melanoma","SC Lung","D-LBCL"};*)


(* ::Input:: *)
(*Legended[PairedBarChart[makeErrors@@#&/@birthList,makeErrors@@#&/@growthList,*)
(*ChartStyle->ColorData[11,"ColorList"],*)
(*FrameLabel->{"A","B"},*)
(*PlotRange->All]*)
(*,Placed[*)
(*SwatchLegend[*)
(*colorList,*)
(*labels,*)
(*LabelStyle->{Black,FontSize->16,FontFamily->"Arial"},*)
(*LegendMarkers->"Bubble",*)
(*LegendLayout->"Row"*)
(*]*)
(*,*)
(*Top]*)
(*]*)


(* ::Input:: *)
(*PairedBarChart[makeErrors@@#&/@birthList,makeErrors@@#&/@growthList,*)
(*ChartStyle->colorList,*)
(*BarSpacing->{5,0,0.5},*)
(*ChartLabels->labels,*)
(*PlotRange->All,*)
(*ImageSize->500,*)
(*ScalingFunctions->"Log",*)
(*FrameStyle->Directive[Thickness[0.0035],Black],*)
(*LabelStyle->{FontSize->18,Black,FontFamily->"Arial"},*)
(*Ticks->{{0.001,0.01,0.1},{}}]*)
(**)


(* ::Input:: *)
(*ColorData[37,"ColorList"]*)


(* ::Input:: *)
(*myleg=SwatchLegend[*)
(*colorList,*)
(*labels,*)
(*LabelStyle->{Black,FontSize->16,FontFamily->"Arial"},*)
(*LegendMarkers->"Bubble",*)
(*LegendLayout->{"Row",3},*)
(*LegendFunction->"Frame"*)
(*]*)


(* ::Input:: *)
(*Export[NotebookDirectory[]<>"../../legend_3Rows.pdf",myleg];*)


(* ::Input:: *)
(*{\[Mu]x,\[Sigma]x} ={0.192, (1.667-0.015)/4};*)
(*{\[Mu]y,\[Sigma]y} ={0.007, (0.033-0.003)/4};*)
(*{\[Mu],\[Sigma]}={Log[\[Mu]x^2/Sqrt[\[Mu]x^2+\[Sigma]x^2]],Log[1+\[Sigma]x^2/\[Mu]x^2]}*)


(* ::Input:: *)
(*Plot[PDF[LogNormalDistribution[-2.514,1.314],x],{x,0,0.5},PlotRange->All]*)


(* ::Input:: *)
(*Mean[LogNormalDistribution[a,b]]*)
(*StandardDeviation[LogNormalDistribution[a,b]]*)


(* ::Input:: *)
(*NSolve[{a+b^2/2==Log[\[Mu]y],(2a+b^2)+Log[-1+E^b^2]==2Log[\[Sigma]y],b>0},{a,b}][[1]]*)


(* ::Input:: *)
(*lognormal[\[Mu]_,min_,max_]:=Module[*)
(*{soln,\[Sigma]},*)
(*\[Sigma]=(max-min)/4;*)
(*soln =Quiet@ NSolve[{a+b^2/2==Log[\[Mu]],(2a+b^2)+Log[-1+E^b^2]==2Log[\[Sigma]],b>0},{a,b}][[1]];*)
(*{a,b}/.soln*)
(*]*)


(* ::Input:: *)
(*birth=lognormal@@#&/@birthList*)
(*netgrowth=lognormal@@#&/@growthList*)


(* ::Input:: *)
(*Plot[Evaluate[PDF[LogNormalDistribution@@#,x]&/@netgrowth],{x,0,0.5},PlotRange->{All,{0,10}}]*)


(* ::Input:: *)
(*drawRatio[bparams_,rparams_]:=Module[*)
(*{},*)
(*Max[0,1-RandomVariate[LogNormalDistribution@@rparams]/RandomVariate[LogNormalDistribution@@bparams]]*)
(*]*)


(* ::Input:: *)
(*randData=Table[drawRatio[birth[[n]],netgrowth[[n]]]&/@Range[200],{n,1,Length[birth]}];*)


(* ::Input:: *)
(*BoxWhiskerChart[randData,"Outliers",*)
(*ChartStyle->colorList,*)
(*BarSpacing->Large,*)
(*ChartLabels->Placed[labels,Axis,Rotate[#,\[Pi]/4]&],*)
(*FrameLabel->{"","\[Mu]/\[Lambda]"},*)
(*FrameStyle->Directive[Thickness[0.0035],Black],*)
(*LabelStyle->{FontSize->14,Black,FontFamily->"Arial"},*)
(*ImageSize->500]*)


(* ::Input:: *)
(*data = Import[NotebookDirectory[]<>"Predicted_LaggingChrPct_TCGA.txt","Data"]*)


(* ::Input:: *)
(*header=data[[1]];*)
(*rawdata=data[[2;;,All]];*)


(* ::Input:: *)
(*Position[rawdata[[All,1]],s_String/;StringMatchQ[s,"brca*",IgnoreCase->True]]//Flatten;*)


(* ::Input:: *)
(*rawdata[[1101,1]]*)


(* ::Input:: *)
(*cancerTypes =DeleteDuplicates[StringTake[rawdata[[All,1]],4]]*)


(* ::Input:: *)
(*locations=Flatten[Position[rawdata[[All,1]],s_String/;StringMatchQ[s,#<>"*",IgnoreCase->True]]]&/@cancerTypes;*)


(* ::Input:: *)
(*BoxWhiskerChart[rawdata[[#,3]]/100&/@locations,"Outliers",*)
(*ChartStyle->colorList,*)
(*BarSpacing->Large,*)
(*ChartLabels->Placed[cancerTypes,Axis,Rotate[#,\[Pi]/4]&],*)
(*FrameLabel->{"","\[Beta]"},*)
(*FrameStyle->Directive[Thickness[0.0035],Black],*)
(*LabelStyle->{FontSize->14,Black,FontFamily->"Arial"},*)
(*ImageSize->500]*)


(* ::Input:: *)
(*labels*)


(* ::Input:: *)
(*DiscretePlot[rawdata[[#,3]]/100&/@locations]*)


(* ::Input:: *)
(*leg2="n_"<>ToString[#]&/@nvals;*)


(* ::Input:: *)
(*(*Export[NotebookDirectory[]<>"criticalcurve_newBetafunc_betaconstant_4Ks.csv",critTable[[#]],"CSV"]&/@Range[1]*)*)
(*Export[NotebookDirectory[]<>"criticalcurve_"<>ToString[leg2[[#]]]<>".csv",critTable[[#]],"CSV"]&/@Range[4]*)


(* ::Input:: *)
(*leg2*)


(* ::Input:: *)
(*critTable2=Table[{\[Delta],critDelta[\[Delta],#]},{\[Delta],1-1*^-3,1-1*^-5,2*^-5}]&/@nvals;*)


(* ::Input:: *)
(*critTable2;*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*Log[critTable2[[All,-1]]]*)


(* ::Input:: *)
(*ListLogLogPlot[critTable2]*)
