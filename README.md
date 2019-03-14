isogeometric_application (IA) (c) Hoang-Giang Bui, 2014, 2015, 2016, 2017, 2018, 2019, Ruhr University Bochum, is an extension of KRATOS Multiphysics framework to assist isogeometric analysis with multipatch NURBS, hierarchical B-Splines and T-Splines. IA uses the Bezier decomposition concept to retain the local characteristics of the finite element. IA is designed to be a complete IGA modeller and employ the Kratos infrastructure to perform analysis.

NURBS:
 + NUBRS patch/multipatch creation. All interface rotation are handled (2D/3D).
 + Support NURBS file format v2.1
 + Utilities to create standard geometry (line, arc, rectangle, ring).
 + Creation of surface by lofting two lines.
 + Creation of volume by lofting two surfaces.
 + Creation of sweep volume.
 + Creation of bending strip patch.
 + h/p/k-refinement. Control values are transferred automatically.
 + Post-processing in GiD (2D/3D, multipatch, using Lagrange mesh).
 + Post-processing in Glvis (single patch).
 + Bezier extraction.

Hierarchical B-Splines:
 + Generation of hierarchical B-Splines patch/multipatch from NURBS patch/multipatch. The conforming of hierarchical B-Splines during refinement is preserved.
 + Local refinement. Control values are transferred automatically.
 + Bezier extraction.

T-Splines:
 + Read-in T-Splines.
 + Bezier extraction.

