:Evaluate: BeginPackage["CHESSNativeBackend`Link`"]
:Evaluate: Begin["`Private`"]

:Begin:
:Function: CHESSNativeLoad
:Pattern: CHESSNativeLoad[file_String]
:Arguments: {file}
:ArgumentTypes: {String}
:ReturnType: Manual
:End:

:Begin:
:Function: CHESSNativeRun
:Pattern: CHESSNativeRun[boundary_String, layers_Integer, columns_Integer, cacheStates_Integer]
:Arguments: {boundary, layers, columns, cacheStates}
:ArgumentTypes: {String, Integer, Integer, Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: CHESSNativePolynomialRun
:Pattern: CHESSNativePolynomialRun[boundary_String, layers_Integer, columns_Integer, cacheStates_Integer]
:Arguments: {boundary, layers, columns, cacheStates}
:ArgumentTypes: {String, Integer, Integer, Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: CHESSNativeFetch
:Pattern: CHESSNativeFetch[mode_Integer, selectedOrder_Integer]
:Arguments: {mode, selectedOrder}
:ArgumentTypes: {Integer, Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: CHESSNativeClear
:Pattern: CHESSNativeClear[]
:Arguments: {}
:ArgumentTypes: {}
:ReturnType: Manual
:End:

:Evaluate: End[]
:Evaluate: EndPackage[]
