function P = defaultRun(varargin)
%This function can be run with no input arguments to run a single triplet.
%The parameters for the triplet are defined within the script. If the
%script is called with any input arguments, it will return a structure P.
%This structure can be input into any of the three FEA solvers and they
%will run that problem. If any function is called without input arguments
%they will call this function and use its output in place of input
%arguments.

global L

%P is an input/output structure. It is capable of storing all of the
%parameters and results for a given problem.
P = struct('numCellsHorz',[],'numCellsVert',[],...%Major problem parameters
           'cellHeight',[],'cellWidth',[],'ligThick',[],...%Major problem parameters
           'latticeType',[],'probType',[], ... %Major problem parameters
           'offset',[],...
           'nelx',[],'nely',[],'plots',[], ...
           'AR',[],... %end of solver parameters
           'strainEnergyCont',[], ... %Raw micropolar results
           'strainEnergyClass',[],... 'UCStrainEnergyClass',[], ...
           'globalStrainEnergyError',[],...
           'globalClassEnergyError',[],...
           'ligLength',[],... %Calculated parameters that describe the lattice.
           'useNewProps',true,... %If true, the continuum solver will use the new correct properties from generalized continuum modeling.
           'MPEffectSize',[]);
P.probType = 'stretchY'; %This sets the category of boundary conditions.
%Common options include:
%'stretchX','stretchY','stretchYfree','transverseY','transverseX','transverseYhinge','curveX','curveY','halfspace',

%%

P.plots = 'displacement';
%P.plots = 'stress,displacement,IEN,blackLattice';
%This line will dictate which plotting options are used as the code runs.
%To include more than one plot option, have a single string where the
%options are separated by commas.
%Options: 'stress,displacement,IEN,blackLattice'
%WARNING: If Matlab opens up hundreds of plots, you can bring your computer
%to a grinding halt.

P.nelx = 8;%NUmber of finite elements in the x and y directions. Default is 35.
P.nely = 8;
P.numCellsHorz = 2; %Number of cells in the horizontal direction.
P.numCellsVert = 2; %Number of cells in the vertical direction.
P.cellHeight = L/P.numCellsVert;
P.cellWidth = L/P.numCellsHorz;



end