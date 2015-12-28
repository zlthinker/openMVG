
ADD_DEFINITIONS(-DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

ADD_EXECUTABLE(openMVG_sample_sfm_test sequential_sfm.cpp)
TARGET_LINK_LIBRARIES(openMVG_sample_sfm_test
        openMVG_image
        openMVG_features
        openMVG_matching
        openMVG_system
        openMVG_sfm
        openMVG_multiview
        openMVG_lInftyComputerVision
        stlplus
        vlsift
        easyexif)

SET_PROPERTY(TARGET openMVG_sample_sfm_test PROPERTY FOLDER OpenMVG/Samples)
