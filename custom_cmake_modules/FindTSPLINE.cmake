SET(TSPLINE_INCLUDE_DIR "${TSPLINE_DIR}/include")

IF(TSPLINE_INCLUDE_DIR)
    SET(TSPLINE_FOUND "YES")
    FIND_LIBRARY(TSPLINE_NEWMAT_LIBRARY newmat "${TSPLINE_DIR}/lib/x86/release" "${TSPLINE_DIR}/lib")
    FIND_LIBRARY(TSPLINE_RHINO_LIBRARY rhino "${TSPLINE_DIR}/lib/x86/release" "${TSPLINE_DIR}/lib")
    FIND_LIBRARY(TSPLINE_LIBRARY tspline "${TSPLINE_DIR}/lib/x86/release" "${TSPLINE_DIR}/lib")
    SET(TSPLINE_LIBRARIES ${TSPLINE_RHINO_LIBRARY} ${TSPLINE_LIBRARY} ${TSPLINE_NEWMAT_LIBRARY})
    message("TSPLINE found")
ELSE()
    message("finding TSPLINE failed, please try to set the var TSPLINE_DIR")
ENDIF()

