//
// Block structured mesh for 2D flow around cylinder
//
//////////////////////////////////////////////////////////////////////
// Parameter section
//////////////////////////////////////////////////////////////////////
// I use open cascade
SetFactory("OpenCASCADE");

// mesh parameters
cyl_rad = DefineNumber[ 0.5 , Min 0.25, Max 1.0, Step 0.05, Name "Cylinder radius" ];
box_square = DefineNumber[ 3.0 , Min 2.0, Max 10.0, Step 1.0, Name "Size of a square region around cylinder" ];
inner_box = DefineNumber[ 1.0 , Min 2.0, Max 10.0, Step 1.0, Name "Size of inner square " ];
box_min_x = DefineNumber[ -15.0 , Min -30.0, Max -10.0, Step 1.0, Name "Domain minimum x" ];
box_max_x = DefineNumber[  35.0 , Min  10.0, Max  50.0, Step 1.0, Name "Domain maximum x" ];
box_width = DefineNumber[ 30.0 , Min 20.0, Max 50.0, Step 1.0, Name "Domain with" ];

// mesh parameters
// this will be the resolution in the fine area
srf_nsplit = DefineNumber[ 8 , Min 2, Max 20, Step 1, Name "Numer of splits; cylinder surface" ];



// number of splits at the cylinder surface from the wake side
wsrf_nsplit = DefineNumber[ 6 , Min 2, Max 20, Step 1, Name "Numer of splits; cylinder surface from the wake side" ];
// number of splits normal to the cylinder surface
nsrf_nsplit = DefineNumber[ 6 , Min 4, Max 20, Step 1, Name "Numer of splits; normal to the cylinder surface" ];
// number of splits in front of cylinder
front_nsplit = DefineNumber[ 6 , Min 2, Max 20, Step 1,   Name "Numer of splits; front of cylinder" ];
// number of splits in wake
wake_nsplit = DefineNumber[ 15 , Min 2, Max 50, Step 1,   Name "Numer of splits; wake" ];
// spanwise number of splits
spnw_nsplit = DefineNumber[ 10 , Min 2, Max 50, Step 1,   Name "Spanwise numer of splits" ];

// progression normal to the cylinder surface
wnprog = DefineNumber[ 1.2 , Min 1, Max 2, Step 0.1,   Name "Progression normal to the cylinder surface" ];
// progression in the wake region
wwprog = DefineNumber[ 1.05 , Min 1, Max 2, Step 0.01,   Name "Progression in the wake region" ];
// progression in the front region
wfprog = DefineNumber[ 1.08 , Min 1, Max 2, Step 0.01,   Name "Progression in the front region" ];

// Element scale at cylinder surface
cs_el_sc = 0.1;

//////////////////////////////////////////////////////////////////////
// Square region around cylinder
//////////////////////////////////////////////////////////////////////
//// Points
//pts_centre = newp;
//Point(newp) = {0.0,0.0,0.0,cs_el_sc};
//
//x_pos = cyl_rad*Sin(Pi/4.0);
//pts_arc_1 = newp;
//Point(newp) = {-x_pos,-x_pos,0.0,cs_el_sc};
//pts_arc_2 = newp;
//Point(newp) = {-x_pos,x_pos,0.0,cs_el_sc};
//
//x_pos = box_square;
//pts_sqr_1 = newp;
//Point(newp) = {-x_pos,-x_pos,0.0,cs_el_sc};
//pts_sqr_2 = newp;
//Point(newp) = {-x_pos,x_pos,0.0,cs_el_sc};
//
//// Curves
//list() = {newl};
//Circle(newl) = {pts_arc_1,pts_centre,pts_arc_2};
//list() += {newl};
//Line(newl) = {pts_arc_2,pts_sqr_2};
//list() += {newl};
//Line(newl) = {pts_sqr_2,pts_sqr_1};
//list() += {newl};
//Line(newl) = {pts_sqr_1,pts_arc_1};
//
//// Surfaces
//surface_list() = {};
//crvl = newll;
//Curve Loop (crvl) = {list()};
//surface_list() += {news};
//Surface(news) = {crvl};
//
//
//// Replicate rotate
//langle = 0.5*Pi;
//surface_list() += Rotate {{0.0,0.0,1.0},{0.0,0.0,0.0},langle} { Duplicata{ Surface{surface_list(0)}; }};
//surface_list() += Rotate {{0.0,0.0,1.0},{0.0,0.0,0.0},-langle} { Duplicata{ Surface{surface_list(0)}; }};
//surface_list() += Rotate {{0.0,0.0,1.0},{0.0,0.0,0.0},2*langle} { Duplicata{ Surface{surface_list(0)}; }};


// this will just be a square for us...
// Points
pts_centre = newp;
Point(newp) = {0.0,0.0,0.0,cs_el_sc};
x_pos = box_square;
pts_sqr_NW = newp;
Point(newp) = {-x_pos,x_pos,0.0,cs_el_sc};
pts_sqr_NE = newp;
Point(newp) = {x_pos,x_pos,0.0,cs_el_sc};
pts_sqr_SE = newp;
Point(newp) = {x_pos,-x_pos,0.0,cs_el_sc};
pts_sqr_SW = newp;
Point(newp) = {-x_pos,-x_pos,0.0,cs_el_sc};

// the guy in the middle
inner = inner_box;
NW = newp;
Point(newp) = {-inner,inner,0.0,cs_el_sc};
NE = newp;
Point(newp) = {inner,inner,0.0,cs_el_sc};
SE = newp;
Point(newp) = {inner,-inner,0.0,cs_el_sc};
SW = newp;
Point(newp) = {-inner,-inner,0.0,cs_el_sc};

// Curves (outer)
list() = {newl};
Line(newl) = {pts_sqr_NW,pts_sqr_NE};
list() += {newl};
Line(newl) = {pts_sqr_NE,pts_sqr_SE};
list() += {newl};
Line(newl) = {pts_sqr_SE,pts_sqr_SW};
list() += {newl};
Line(newl) = {pts_sqr_SW,pts_sqr_NW};

// Curves (inner)
list_in() = {newl};
Line(newl) = {NW,NE};
list_in() += {newl};
Line(newl) = {NE,SE};
list_in() += {newl};
Line(newl) = {SE,SW};
list_in() += {newl};
Line(newl) = {SW,NW};


// Curves (between)
list_in() += {newl};
Line(newl) = {pts_sqr_NW,NW};
list_in() += {newl};
Line(newl) = {pts_sqr_NE,NE};
list_in() += {newl};
Line(newl) = {pts_sqr_SE,SE};
list_in() += {newl};
Line(newl) = {pts_sqr_SW,SW};

// Surfaces
surface_list() = {};

// little inner guy
crvl = newll;
Curve Loop (crvl) = {list_in(0),list_in(1),list_in(2),list_in(3)};
surface_list() += {news};
Surface(news) = {crvl};


// top
crvl = newll;
Curve Loop (crvl) = {list_in(0),list_in(4),list(0),list_in(5)};
surface_list() += {news};
Surface(news) = {crvl};

// right
crvl = newll;
Curve Loop (crvl) = {list_in(1),list_in(5),list(1),list_in(6)};
surface_list() += {news};
Surface(news) = {crvl};

// bottom
crvl = newll;
Curve Loop (crvl) = {list_in(2),list_in(6),list(2),list_in(7)};
surface_list() += {news};
Surface(news) = {crvl};


// left
crvl = newll;
Curve Loop (crvl) = {list_in(3),list_in(7),list(3),list_in(4)};
surface_list() += {news};
Surface(news) = {crvl};

// Remove duplicates
Coherence;

// Extract edges
list() = Unique(Abs(Boundary { Surface{surface_list()}; }));
//For il In {0:# list() - 1}
//   Printf("list %g %g", il,list(il));
//EndFor

// Extrude
// -x
ltmp[] = Extrude{box_square+box_min_x,0,0} { Curve{list(3)}; };
surface_list() += {ltmp[1]};
// +x
ltmp[] = Extrude{box_max_x-box_square,0,0} { Curve{list(1)}; };
surface_list() += {ltmp[1]};
// -y
ltmp[] = Extrude{0,box_square-box_width*0.5,0,0} { Curve{list(2)}; };
surface_list() += {ltmp[1]};
// +y
ltmp[] = Extrude{0,box_width*0.5-box_square,0,0} { Curve{list(0)}; };
surface_list() += {ltmp[1]};

// Remove duplicates
Coherence;

// Extract edges
list() = Unique(Abs(Boundary { Surface{surface_list()}; }));
For il In {0:# list() - 1}
   Printf("list HARRY %g %g", il,list(il));
EndFor

// Extrude
// -x -y
ltmp[] = Extrude{0,box_square-box_width*0.5,0,0} { Curve{13}; };
surface_list() += {ltmp[1]};
// -x +y
ltmp[] = Extrude{0,box_width*0.5-box_square,0,0} { Curve{14}; };
surface_list() += {ltmp[1]};
// +x -y
ltmp[] = Extrude{0,box_square-box_width*0.5,0,0} { Curve{17}; };
surface_list() += {ltmp[1]};
// +x +y
ltmp[] = Extrude{0,box_width*0.5-box_square,0,0} { Curve{16}; };
surface_list() += {ltmp[1]};

// Remove duplicates
Coherence;

// Extract edges
list() = Unique(Abs(Boundary { Surface{surface_list()}; }));
For il In {0:# list() - 1}
   Printf("list %g %g", il,list(il));
EndFor

//////////////////////////////////////////////////////////////////////
// Physical properties section
//////////////////////////////////////////////////////////////////////
// Volume
Physical Surface("Fluid",1) = {surface_list()};

// Wall
// Physical Curve("Wal",2) = {list(3),list(6),list(9),list(11)};

// Inflow
//Physical Curve("v  ",1) = {list(22),list(6),list(20)};
Physical Curve("v  ",1) = {27,15,25};

// Outflow
// Physical Curve("o  ",2) = {list(18),list(9),list(16)};
Physical Curve("o  ",2) = {31,18,29};

// Top periodic
//Physical Curve("Pet",3) = {list(23),list(15),list(19)};
Physical Curve("Pet",3) = {28,24,32};

// Bottom periodic
//Physical Curve("Peb",4) = {list(21),list(12),list(17)};
Physical Curve("Peb",4) = {26,21,30};

//////////////////////////////////////////////////////////////////////
// Transfinite division section
//////////////////////////////////////////////////////////////////////
// Edge division
//
// Transfinite Curve{list(3),list(6),list(9),list(1),list(5),list(7),list(14),list(20),list(23)} = (srf_nsplit+1) Using Progression 1;

// Transfinite Curve{list(11),list(10),list(17)} = (wsrf_nsplit+1) Using Progression 1;

// this is all the inside stuff
// bit of a judgement call here...
// Harry
// lets say 2x the surface splitting??
//Transfinite Curve{list(6),list(3),list(1),list(9),list(12),list(2),list(0),list(15)} = (2*srf_nsplit+1) Using Progression 1;
//Transfinite Curve{14,3,7,5,1,17,23,0,4,6,2,20} = (2*srf_nsplit+1) Using Progression 1;
Transfinite Curve{15,4,8,6,2,18,24,1,5,7,3,21} = (srf_nsplit) Using Progression 1;

// this is the outer mesh
//Transfinite Curve{list(22),list(13),list(14),list(18),list(20),list(11),list(10),list(16)} = (spnw_nsplit+1) Using Progression 1;
Transfinite Curve{27,22,23,31,29,19,20,25} = (spnw_nsplit+1) Using Progression 1;

//Transfinite Curve{list(21),list(4),list(5),list(23)} = (front_nsplit+1) Using Progression wfprog;
Transfinite Curve{28,14,13,26} = (front_nsplit+1) Using Progression wfprog;

//Transfinite Curve{list(17),list(8),list(7),list(19)} = (wake_nsplit+1) Using Progression wwprog;
Transfinite Curve{32,16,17,30} = (wake_nsplit+1) Using Progression wwprog;

// no cylinder
// Transfinite Curve{-list(0),list(2),-list(4),list(8)} = (nsrf_nsplit+1) Using Progression wnprog;
Transfinite Curve{-9,-10,-11,-12} = (nsrf_nsplit+1) Using Progression wnprog;

// Surface division
Transfinite Surface{surface_list()};
Recombine Surface{surface_list()};


//////////////////////////////////////////////////////////////////////
// Meshing section
//////////////////////////////////////////////////////////////////////
Mesh 1;
Mesh 2;
Mesh 3;

SetOrder 2;

RenumberMeshElements;

//////////////////////////////////////////////////////////////////////
// Mesh saving section
//////////////////////////////////////////////////////////////////////
Mesh.Format = 1;
Mesh.MshFileVersion = 2.2;
Mesh.SaveAll = 0;
Mesh.Binary = 0;

Save "brink_cyl.msh";

