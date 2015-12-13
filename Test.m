

%Over one cycle
%Each measurement (observation) will create the following structures
struct_AngularPosition = struct('Longitude',0, 'Latitude',0);%Coordinates of earth
struct_GridPosition = struct('Row',0,'Column',0);%Coordinates on grid (Row-Column format)
struct_SSH = struct('SSH', 0);%SSH data
%Structure for single measurement
struct_ToGUIComponent = struct('EarthLocation', struct_AngularPosition,'GridPosition',struct_GridPosition,'SSH',struct_SSH);
%Allocating memory for multiple measurements
n_Observations = 100;
allocate_Struct = repmat(struct('EarthLocation',struct('Longitude',0,'Latitude',0),'GridPosition',struct('Row',0','Column',0),'SSH',struct('SSH',0)),n_Observations,1);

struct_CurrentAngularPosition = allocate_Struct(1).EarthLocation;
struct_CurrentAngularPosition.Longitude = 1
struct_CurrentAngularPosition.Latitude = 1
struct_SSH = allocate_Struct(1).SSH;