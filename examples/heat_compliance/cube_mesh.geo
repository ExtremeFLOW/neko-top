SetFactory("OpenCASCADE");

// ------------------------------------------------------------
// Parameters
// ------------------------------------------------------------
L  = 3;        // Cube side length
N  = 11;        // Default subdivisions per coarse cell (override with -setnumber N)

// Cell size so that everything is cubes
h  = L/(3*N);  // => 3N cells per cube edge
dx = L/3;
lc = h;
Nz = 3*N;      // layers in z so dz = h

// Tell Gmsh to use 2nd order, *full* elements (quad9 / hex27)
Mesh.ElementOrder           = 2;
Mesh.SecondOrderIncomplete  = 0;  // 0 => full order (9-node quads, 27-node hexes)
Mesh.RecombineAll           = 1;  // hexes instead of tets
// Mesh.HighOrderOptimize   = 1;  // optional smoothing

// ------------------------------------------------------------
// Bottom grid: points (i=0..3, j=0..3) at z = 0
// Indexing: idx = i + 4*j
// ------------------------------------------------------------
For j In {0:3}
  For i In {0:3}
    idx = i + 4*j;
    p[idx] = newp;
    Point(p[idx]) = { i*dx, j*dx, 0, lc };
  EndFor
EndFor

// ------------------------------------------------------------
// Lines along x (horizontal) and y (vertical) on bottom
// ------------------------------------------------------------

// Horizontal lines: hx[]
k = 0;
For j In {0:3}
  For i In {0:2}
    idL = i     + 4*j;
    idR = (i+1) + 4*j;
    hx[k] = newl;
    Line(hx[k]) = { p[idL], p[idR] };
    k += 1;
  EndFor
EndFor

// Vertical lines: hy[]
k = 0;
For i In {0:3}
  For j In {0:2}
    idB = i + 4*j;
    idT = i + 4*(j+1);
    hy[k] = newl;
    Line(hy[k]) = { p[idB], p[idT] };
    k += 1;
  EndFor
EndFor

// Subdivide each coarse segment into N elements (N+1 points)
For k In {0:#hx[]-1}
  Transfinite Line{ hx[k] } = N + 1;
EndFor

For k In {0:#hy[]-1}
  Transfinite Line{ hy[k] } = N + 1;
EndFor

// ------------------------------------------------------------
// Build 3x3 bottom patches and extrude each
// ------------------------------------------------------------
bottom_other[] = {};
top_surfs[]    = {};
side_surfs[]   = {};
vols[]         = {};
center_surf    = 0;

For j In {0:2}
  For i In {0:2}
    // Map patch (i,j) to its four boundary lines:
    // hx index: row j, segment i       -> i + 3*j
    // hy index: column i, segment j    -> j + 3*i
    south = hx[i + 3*j];
    east  = hy[j + 3*(i+1)];
    north = hx[i + 3*(j+1)];
    west  = hy[j + 3*i];

    ll = newll;
    Line Loop(ll) = { south, east, -north, -west };

    s = news;
    Plane Surface(s) = { ll };
    Transfinite Surface{s};
    Recombine Surface{s};

    // Central patch (middle of 3x3) -> BC id 2
    If ( (i == 1) && (j == 1) )
      center_surf = s;
    Else
      bottom_other[] += { s };
    EndIf

    // Extrude to build full cube. Nz layers => cubic cells.
    out[] = Extrude {0, 0, L} {
      Surface{s};
      Layers{ Nz };
      Recombine;
    };

    // out[0] : top surface of this column
    // out[1] : volume
    // out[2..5] : side surfaces (south, east, north, west)
    top_surfs[] += { out[0] };
    vols[]      += { out[1] };

    // Collect only external side surfaces
    If (j == 0)
      side_surfs[] += { out[2] }; // y = 0 (south)
    EndIf
    If (i == 2)
      side_surfs[] += { out[3] }; // x = 3 (east)
    EndIf
    If (j == 2)
      side_surfs[] += { out[4] }; // y = 3 (north)
    EndIf
    If (i == 0)
      side_surfs[] += { out[5] }; // x = 0 (west)
    EndIf
  EndFor
EndFor

// ------------------------------------------------------------
// Transfinite volumes (structured hex grid)
// ------------------------------------------------------------
Transfinite Volume{ vols[] };

// ------------------------------------------------------------
// Physical groups *with names* (important for gmsh2nek)
// ------------------------------------------------------------

// Volume
Physical Volume("fluid", 1) = { vols[] };

// Bottom center patch: BC id 2
Physical Surface("bottom_center", 2) = { center_surf };

// Everything else on external boundary: BC id 1
all_bc1[] = {};
all_bc1[] += { bottom_other[] };
all_bc1[] += { top_surfs[] };
all_bc1[] += { side_surfs[] };

Physical Surface("outer_walls", 1) = { all_bc1[] };


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

Save "cube_mesh.msh";
