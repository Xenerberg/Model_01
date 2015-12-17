%---------------------------------------------------------------------------
%  MATLAB function for configuration set generated on 17-Dec-2015 12:31:39
%  MATLAB version: 8.3.0.532 (R2014a)
%---------------------------------------------------------------------------

function cs=config_test()

cs = Simulink.ConfigSet;
 
% Original configuration set version: 1.14.2
if cs.versionCompare('1.14.2') < 0
    error('Simulink:MFileVersionViolation', 'The version of the target configuration set is older than the original configuration set.');
end
 
% Original environment character encoding: windows-1252
if ~strcmpi(get_param(0, 'CharacterEncoding'), 'windows-1252')
    warning('Simulink:EncodingUnMatched', 'The target character encoding (%s) is different from the original (%s).',  get_param(0, 'CharacterEncoding'), 'windows-1252');
end
 
% Do not change the order of the following commands. There are dependencies between the parameters.
hdlcoderui_hdlcc_Successful=true;
try
   componentCC = hdlcoderui.hdlcc;
   cs.attachComponent(componentCC);
catch ME
   warning('Simulink:ConfigSet:AttachComponentError', ME.message);
   hdlcoderui_hdlcc_Successful=false;
end

cs.set_param('SystemTargetFile', 'grt.tlc');   % System target file
cs.set_param('GenerateComments', 'on');   % Include comments 
cs.set_param('SimulinkBlockComments', 'on');   % Simulink block / Stateflow object comments 
cs.set_param('ShowEliminatedStatement', 'off');   % Show eliminated blocks 
cs.set_param('ForceParamTrailComments', 'off');   % Verbose comments for SimulinkGlobal storage class 
cs.set_param('GenerateReport', 'off');   % Create code generation report 
cs.set_param('LaunchReport', 'off');   % Open report automatically 
cs.set_param('InlineParams', 'off');   % Inline parameters 
cs.set_param('OptimizeBlockIOStorage', 'on');   % Signal storage reuse 
cs.set_param('ExpressionFolding', 'on');   % Eliminate superfluous local variables (expression folding) 
cs.set_param('UpdateModelReferenceTargets', 'IfOutOfDateOrStructuralChange');   % Rebuild 
cs.set_param('ProdHWDeviceType', '32-bit Generic');   % Device type 
cs.set_param('EnableMemcpy', 'on');   % Use memcpy for vector assignment 
cs.set_param('SaveFinalState', 'off');   % Final states 
cs.set_param('SignalLogging', 'on');   % Signal logging 
cs.set_param('TargetUnknown', 'off');   % Target unknown 
cs.set_param('DSMLogging', 'on');   % Data stores 
cs.set_param('EnableConcurrentExecution', 'off');   % EnableConcurrentExecution 
cs.set_param('LimitDataPoints', 'on');   % Limit data points to last 
cs.set_param('LoadExternalInput', 'off');   % Input 
cs.set_param('LoadInitialState', 'off');   % Initial state 
cs.set_param('ReturnWorkspaceOutputs', 'off');   % Save simulation output as single object 
cs.set_param('SaveOutput', 'on');   % Output 
cs.set_param('SaveState', 'off');   % States 
cs.set_param('SaveTime', 'on');   % Time 
cs.set_param('SupportNonFinite', 'on');   % Support non-finite numbers 
cs.set_param('MatFileLogging', 'on');   % MAT-file logging 
cs.set_param('GRTInterface', 'off');   % Classic call interface 
cs.set_param('ModelReferenceMinAlgLoopOccurrences', 'off');   % Minimize algebraic loop occurrences 
cs.set_param('ProdEqTarget', 'on');   % Test hardware is the same as production hardware 
cs.set_param('SolverType', 'Variable-step');   % Type 
cs.set_param('AutoInsertRateTranBlk', 'off');   % Automatically handle rate transition for data transfer 
cs.set_param('ExtMode', 'off');   % Interface 
 
% The following commands do not have dependencies.
cs.set_param('Description', '');  % Description
cs.set_param('Name', 'Configuration');  % Name

% Solver pane
cs.set_param('AbsTol', 'auto');   % Absolute tolerance 
cs.set_param('ConsecutiveZCsStepRelTol', '10*128*eps');   % Time tolerance 
cs.set_param('InitialStep', 'auto');   % Initial step size 
cs.set_param('MaxConsecutiveMinStep', '1');   % Number of consecutive min steps 
cs.set_param('MaxConsecutiveZCs', '1000');   % Number of consecutive zero crossings 
cs.set_param('MaxStep', 'auto');   % Max step size 
cs.set_param('MinStep', 'auto');   % Min step size 
cs.set_param('PositivePriorityOrder', 'off');   % Higher priority value indicates higher task priority 
cs.set_param('RelTol', '1e-3');   % Relative tolerance 
cs.set_param('ShapePreserveControl', 'DisableAll');   % Shape preservation 
cs.set_param('Solver', 'ode45');   % Solver 
cs.set_param('StartTime', '0.0');   % Start time 
cs.set_param('StopTime', '10.0');   % Stop time 
cs.set_param('ZcThreshold', 'auto');   % Signal threshold 
cs.set_param('ZeroCrossAlgorithm', 'Nonadaptive');   % Algorithm 
cs.set_param('ZeroCrossControl', 'UseLocalSettings');   % Zero-crossing control 

% Data Import/Export pane
cs.set_param('DSMLoggingName', 'dsmout');   % Data stores logging name 
cs.set_param('Decimation', '1');   % Decimation 
cs.set_param('InspectSignalLogs', 'off');   % Record and inspect simulation output 
cs.set_param('MaxDataPoints', '1000');   % Limit data points to last 
cs.set_param('OutputOption', 'RefineOutputTimes');   % Output options 
cs.set_param('OutputSaveName', 'yout');   % Output 
cs.set_param('Refine', '1');   % Refine factor 
cs.set_param('SaveFormat', 'Array');   % Format 
cs.set_param('SignalLoggingName', 'logsout');   % Signal logging name 
cs.set_param('SignalLoggingSaveFormat', 'Dataset');   % Signal logging format 
cs.set_param('TimeSaveName', 'tout');   % Time 

% Optimization pane
cs.set_param('AccelVerboseBuild', 'off');   % Verbose accelerator builds 
cs.set_param('BlockReduction', 'on');   % Block reduction 
cs.set_param('BooleanDataType', 'on');   % Implement logic signals as Boolean data (vs. double) 
cs.set_param('ConditionallyExecuteInputs', 'on');   % Conditional input branch execution 
cs.set_param('DefaultUnderspecifiedDataType', 'double');   % Default for underspecified data type 
cs.set_param('EfficientFloat2IntCast', 'off');   % Remove code from floating-point to integer conversions that wraps out-of-range values 
cs.set_param('EfficientMapNaN2IntZero', 'on');   % Remove code from floating-point to integer conversions with saturation that maps NaN to zero 
cs.set_param('InitFltsAndDblsToZero', 'off');   % Use memset to initialize floats and doubles to 0.0 
cs.set_param('LifeSpan', 'inf');   % Application lifespan (days) 
cs.set_param('SimCompilerOptimization', 'Off');   % Compiler optimization level 
cs.set_param('UseFloatMulNetSlope', 'off');   % Use floating-point multiplication to handle net slope corrections 
cs.set_param('UseIntDivNetSlope', 'off');   % Use integer division to handle net slopes that are reciprocals of integers 

% Optimization:Signals and Parameters pane
cs.set_param('BufferReuse', 'on');   % Reuse local block outputs 
cs.set_param('LocalBlockOutputs', 'on');   % Enable local block outputs 
cs.set_param('MaxStackSize', 'Inherit from target');   % Maximum stack size (bytes) 
cs.set_param('MemcpyThreshold', 64);   % Memcpy threshold (bytes) 
cs.set_param('RollThreshold', 5);   % Loop unrolling threshold 

% Optimization:Stateflow pane
cs.set_param('DataBitsets', 'off');   % Use bitsets for storing Boolean data 
cs.set_param('StateBitsets', 'off');   % Use bitsets for storing state configuration 

% Diagnostics Solver pane
cs.set_param('AlgebraicLoopMsg', 'warning');   % Algebraic loop 
cs.set_param('ArtificialAlgebraicLoopMsg', 'warning');   % Minimize algebraic loop 
cs.set_param('BlockPriorityViolationMsg', 'warning');   % Block priority violation 
cs.set_param('ConsistencyChecking', 'none');   % Solver data inconsistency 
cs.set_param('MaxConsecutiveZCsMsg', 'error');   % Consecutive zero crossings violation 
cs.set_param('MinStepSizeMsg', 'warning');   % Min step size violation 
cs.set_param('ModelReferenceExtraNoncontSigs', 'error');   % Extraneous discrete derivative signals 
cs.set_param('SimStateInterfaceChecksumMismatchMsg', 'warning');   % SimState interface checksum mismatch 
cs.set_param('SimStateOlderReleaseMsg', 'error');   % SimState object from earlier release 
cs.set_param('SolverPrmCheckMsg', 'warning');   % Automatic solver parameter selection 
cs.set_param('StateNameClashWarn', 'none');   % State name clash 
cs.set_param('TimeAdjustmentMsg', 'none');   % Sample hit time adjusting 
cs.set_param('UnknownTsInhSupMsg', 'warning');   % Unspecified inheritability of sample time 

% Diagnostics:Sample Time pane
cs.set_param('DiscreteInheritContinuousMsg', 'warning');   % Discrete used as continuous 
cs.set_param('InheritedTsInSrcMsg', 'warning');   % Source block specifies -1 sample time 
cs.set_param('MultiTaskCondExecSysMsg', 'error');   % Multitask conditionally executed subsystem 
cs.set_param('MultiTaskRateTransMsg', 'error');   % Multitask rate transition 
cs.set_param('SigSpecEnsureSampleTimeMsg', 'warning');   % Enforce sample times specified by Signal Specification blocks 
cs.set_param('SingleTaskRateTransMsg', 'none');   % Single task rate transition 
cs.set_param('TasksWithSamePriorityMsg', 'warning');   % Tasks with equal priority 

% Diagnostics:Data Validity pane
cs.set_param('ArrayBoundsChecking', 'none');   % Array bounds exceeded 
cs.set_param('AssertControl', 'UseLocalSettings');   % Model Verification block enabling 
cs.set_param('CheckMatrixSingularityMsg', 'none');   % Division by singular matrix 
cs.set_param('IntegerOverflowMsg', 'warning');   % Detect overflow 
cs.set_param('MultiTaskDSMMsg', 'error');   % Multitask data store 
cs.set_param('ParameterDowncastMsg', 'error');   % Detect downcast 
cs.set_param('ParameterOverflowMsg', 'error');   % Detect overflow 
cs.set_param('ParameterPrecisionLossMsg', 'warning');   % Detect precision loss 
cs.set_param('ParameterTunabilityLossMsg', 'warning');   % Detect loss of tunability 
cs.set_param('ParameterUnderflowMsg', 'none');   % Detect underflow 
cs.set_param('RTPrefix', 'error');   % "rt" prefix for identifiers 
cs.set_param('ReadBeforeWriteMsg', 'UseLocalSettings');   % Detect read before write 
cs.set_param('SignalInfNanChecking', 'none');   % Inf or NaN block output 
cs.set_param('SignalRangeChecking', 'none');   % Simulation range checking 
cs.set_param('SignalResolutionControl', 'UseLocalSettings');   % Signal resolution 
cs.set_param('UnderSpecifiedDataTypeMsg', 'none');   % Underspecified data types 
cs.set_param('UnderspecifiedInitializationDetection', 'Simplified');   % Underspecified initialization detection 
cs.set_param('UniqueDataStoreMsg', 'none');   % Duplicate data store names 
cs.set_param('WriteAfterReadMsg', 'UseLocalSettings');   % Detect write after read 
cs.set_param('WriteAfterWriteMsg', 'UseLocalSettings');   % Detect write after write 

% Diagnostics:Type Conversion pane
cs.set_param('FixptConstOverflowMsg', 'none');   % Detect overflow 
cs.set_param('FixptConstPrecisionLossMsg', 'none');   % Detect precision loss 
cs.set_param('FixptConstUnderflowMsg', 'none');   % Detect underflow 
cs.set_param('Int32ToFloatConvMsg', 'warning');   % 32-bit integer to single precision float conversion 
cs.set_param('UnnecessaryDatatypeConvMsg', 'none');   % Unnecessary type conversions 
cs.set_param('VectorMatrixConversionMsg', 'none');   % Vector/matrix block input conversion 

% Diagnostics:Connectivity pane
cs.set_param('BusNameAdapt', 'WarnAndRepair');   % Repair bus selections 
cs.set_param('BusObjectLabelMismatch', 'warning');   % Element name mismatch 
cs.set_param('FcnCallInpInsideContextMsg', 'EnableAllAsError');   % Context-dependent inputs 
cs.set_param('InvalidFcnCallConnMsg', 'error');   % Invalid function-call connection 
cs.set_param('NonBusSignalsTreatedAsBus', 'none');   % Non-bus signals treated as bus signals 
cs.set_param('RootOutportRequireBusObject', 'warning');   % Unspecified bus object at root Outport block 
cs.set_param('SignalLabelMismatchMsg', 'none');   % Signal label mismatch 
cs.set_param('StrictBusMsg', 'ErrorLevel1');   % Bus signal treated as vector 
cs.set_param('UnconnectedInputMsg', 'warning');   % Unconnected block input ports 
cs.set_param('UnconnectedLineMsg', 'warning');   % Unconnected line 
cs.set_param('UnconnectedOutputMsg', 'warning');   % Unconnected block output ports 

% Diagnostics:Compatibility pane
cs.set_param('FrameProcessingCompatibilityMsg', 'warning');   % Block behavior depends on frame status of signal 
cs.set_param('SFcnCompatibilityMsg', 'none');   % S-function upgrades needed 

% Diagnostics:Model Referencing pane
cs.set_param('ModelReferenceCSMismatchMessage', 'none');   % Model configuration mismatch 
cs.set_param('ModelReferenceDataLoggingMessage', 'warning');   % Unsupported data logging 
cs.set_param('ModelReferenceIOMismatchMessage', 'none');   % Port and parameter mismatch 
cs.set_param('ModelReferenceIOMsg', 'none');   % Invalid root Inport/Outport block connection 
cs.set_param('ModelReferenceVersionMismatchMessage', 'none');   % Model block version mismatch 

% Diagnostics:Saving pane
cs.set_param('SaveWithDisabledLinksMsg', 'warning');   % Block diagram contains disabled library links 
cs.set_param('SaveWithParameterizedLinksMsg', 'warning');   % Block diagram contains parameterized library links 

% Diagnostics:Stateflow pane
cs.set_param('SFInvalidInputDataAccessInChartInitDiag', 'warning');   % Invalid input data access in chart initialization 
cs.set_param('SFNoUnconditionalDefaultTransitionDiag', 'warning');   % No unconditional default transitions 
cs.set_param('SFTransitionActionBeforeConditionDiag', 'warning');   % Transition action specified before condition action 
cs.set_param('SFTransitionOutsideNaturalParentDiag', 'warning');   % Transition outside natural parent 
cs.set_param('SFUnconditionalTransitionShadowingDiag', 'warning');   % Transition shadowing 
cs.set_param('SFUndirectedBroadcastEventsDiag', 'warning');   % Undirected event broadcasts 
cs.set_param('SFUnexpectedBacktrackingDiag', 'warning');   % Unexpected backtracking 
cs.set_param('SFUnusedDataAndEventsDiag', 'warning');   % Unused data and events 

% Hardware Implementation pane
cs.set_param('ProdIntDivRoundTo', 'Undefined');   % Signed integer division rounds to 
cs.set_param('ProdLargestAtomicFloat', 'None');   % floating-point 
cs.set_param('ProdLargestAtomicInteger', 'Char');   % integer 

% Model Referencing pane
cs.set_param('EnableParallelModelReferenceBuilds', 'off');   % Enable parallel model reference builds 
cs.set_param('ModelDependencies', '');   % Model dependencies 
cs.set_param('ModelReferenceNumInstancesAllowed', 'Multi');   % Total number of instances allowed per top model 
cs.set_param('ModelReferencePassRootInputsByReference', 'on');   % Pass fixed-size scalar root inputs by value for code generation 
cs.set_param('ParallelModelReferenceMATLABWorkerInit', 'None');   % MATLAB worker initialization for builds 
cs.set_param('PropagateSignalLabelsOutOfModel', 'off');   % Propagate all signal labels out of the model 
cs.set_param('PropagateVarSize', 'Infer from blocks in model');   % Propagate sizes of variable-size signals 

% Simulation Target General pane
cs.set_param('SFSimEcho', 'on');   % Echo expressions without semicolons 
cs.set_param('SFSimEnableDebug', 'on');   % Enable debugging/animation 
cs.set_param('SFSimOverflowDetection', 'on');   % Enable overflow detection (with debugging) 
cs.set_param('SimBuildMode', 'sf_incremental_build');   % Simulation target build mode 
cs.set_param('SimCtrlC', 'on');   % Ensure responsiveness 
cs.set_param('SimGenImportedTypeDefs', 'off');   % Generate typedefs for imported bus and enumeration types 
cs.set_param('SimIntegrity', 'on');   % Ensure memory integrity 

% Simulation Target:Custom Code pane
cs.set_param('SimCustomHeaderCode', '');   % Header file 
cs.set_param('SimCustomInitializer', '');   % Initialize function 
cs.set_param('SimCustomSourceCode', '');   % Source file 
cs.set_param('SimCustomTerminator', '');   % Terminate function 
cs.set_param('SimParseCustomCode', 'on');   % Parse custom code symbols 
cs.set_param('SimUserIncludeDirs', '');   % Include directories 
cs.set_param('SimUserLibraries', '');   % Libraries 
cs.set_param('SimUserSources', '');   % Source files 

% Simulation Target:Symbols pane
cs.set_param('SimReservedNameArray', '');   % Reserved names 

% Code Generation General pane
cs.set_param('BuildConfiguration', 'Faster Builds');    
cs.set_param('CheckMdlBeforeBuild', 'Off');   % Check model before generating code 
cs.set_param('GenCodeOnly', 'off');   % Generate code only 
cs.set_param('ObjectivePriorities', '');   % Prioritized objectives 
cs.set_param('PackageGeneratedCodeAndArtifacts', 'off');   % Package code and artifacts 
cs.set_param('PackageName', '');   % Zip file name 
cs.set_param('TargetLang', 'C');   % Language 
cs.set_param('Toolchain', 'Automatically locate an installed toolchain');   % Toolchain 

% Code Generation:Comments pane
cs.set_param('MATLABSourceComments', 'off');   % MATLAB source code as comments 

% Code Generation:Symbols pane
cs.set_param('MaxIdLength', 31);   % Maximum identifier length 
cs.set_param('ReservedNameArray', '');   % Reserved names 
cs.set_param('UseSimReservedNames', 'off');   % Use the same reserved names as Simulation Target 

% Code Generation:Custom Code pane
cs.set_param('CustomHeaderCode', '');   % Header file 
cs.set_param('CustomInclude', '');   % Include directories 
cs.set_param('CustomInitializer', '');   % Initialize function 
cs.set_param('CustomLibrary', '');   % Libraries 
cs.set_param('CustomSource', '');   % Source files 
cs.set_param('CustomSourceCode', '');   % Source file 
cs.set_param('CustomTerminator', '');   % Terminate function 
cs.set_param('RTWUseSimCustomCode', 'off');   % Use the same custom code settings as Simulation Target 

% Code Generation:Debug pane
cs.set_param('ProfileTLC', 'off');   % Profile TLC 
cs.set_param('RTWVerbose', 'on');   % Verbose build 
cs.set_param('RetainRTWFile', 'off');   % Retain .rtw file 
cs.set_param('TLCAssert', 'off');   % Enable TLC assertion 
cs.set_param('TLCCoverage', 'off');   % Start TLC coverage when generating code 
cs.set_param('TLCDebug', 'off');   % Start TLC debugger when generating code 

% Code Generation:Interface pane
cs.set_param('CodeInterfacePackaging', 'Nonreusable function');   % Code interface packaging 
cs.set_param('CodeReplacementLibrary', 'None');   % Code replacement library 
cs.set_param('GenerateASAP2', 'off');   % Interface 
cs.set_param('LogVarNameModifier', 'rt_');   % MAT-file variable name modifier 
cs.set_param('TargetLangStandard', 'C89/C90 (ANSI)');   % Standard math library 
cs.set_param('UtilityFuncGeneration', 'Auto');   % Shared code placement 

