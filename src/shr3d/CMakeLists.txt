SET(SHR3D_HEADER_FILES
        shr3d.h)

SET(SHR3D_SOURCE_FILES
        shr3d.cpp)

ADD_LIBRARY(SHR3D_LIB STATIC ${SHR3D_HEADER_FILES} ${SHR3D_SOURCE_FILES})
TARGET_LINK_LIBRARIES(SHR3D_LIB
        PUBLIC
        PUBGEO_LIB)

### Add some executables to play with.
### Main exe will showcase all features
ADD_EXECUTABLE(shr3d main.cpp)
TARGET_LINK_LIBRARIES(shr3d
        PUBLIC
        SHR3D_LIB)

IF(TESTING)
    ## This executable will be used to test IO
    ADD_EXECUTABLE(shr3dTest test-main.cpp)
    TARGET_LINK_LIBRARIES(shr3dTest
            PUBLIC
            SHR3D_LIB)
ENDIF()
