:Evaluate: BeginPackage["CHESSNativeBackend`Link`"]
:Evaluate: Begin["`Private`"]

:Begin:
:Function: chess_native_load
:Pattern: ChessNativeLoad[file_String]
:Arguments: {file}
:ArgumentTypes: {String}
:ReturnType: Manual
:End:

:Begin:
:Function: chess_native_run
:Pattern: ChessNativeRun[boundary_String, layers_Integer, columns_Integer, cacheStates_Integer]
:Arguments: {boundary, layers, columns, cacheStates}
:ArgumentTypes: {String, Integer, Integer, Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: chess_native_fetch
:Pattern: ChessNativeFetch[mode_Integer, selectedOrder_Integer]
:Arguments: {mode, selectedOrder}
:ArgumentTypes: {Integer, Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: chess_native_clear
:Pattern: ChessNativeClear[]
:Arguments: {}
:ArgumentTypes: {}
:ReturnType: Manual
:End:

:Evaluate: End[]
:Evaluate: EndPackage[]
