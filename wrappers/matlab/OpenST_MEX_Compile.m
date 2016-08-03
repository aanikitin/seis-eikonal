function OpenST_MEX_Compile()
    mex -largeArrayDims COMPFLAGS='$COMPFLAGS /I "../../include"' ...
        OpenST_LSM3D_mex.c OpenST_MEX.c ...
        ../../lib/Release/openst-0.2-static.lib
    mex -largeArrayDims COMPFLAGS='$COMPFLAGS /I "../../include"' ...
        OpenST_BRT3D_mex.c OpenST_MEX.c ...
        ../../lib/Release/openst-0.2-static.lib
end
