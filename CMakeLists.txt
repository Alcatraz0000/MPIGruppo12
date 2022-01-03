cmake_minimum_required(VERSION 3.9)
project (Contest-MPI
    DESCRIPTION "RadixSort implementation with MPI"
	VERSION 0.1
	LANGUAGES C
)

enable_testing()

SET (INCLUDE_LIBRARY_C RadixSort.c CACHE STRING "Library Name C")
SET (INCLUDE_LIBRARY_H RadixSort.h CACHE STRING "Library Name H")

SET (INCLUDE_LIBRARY_SEQ_C RadixSortSeq.c CACHE STRING "Library Seq Name C")
SET (INCLUDE_LIBRARY_SEQ_H RadixSortSeq.h CACHE STRING "Library Seq Name H")



#create Executable. The OX rappresents the opt option.
add_executable(${PROJECT_NAME}_O1 src/main.c src/${INCLUDE_LIBRARY_C} include/${INCLUDE_LIBRARY_H})
add_executable(${PROJECT_NAME}_O2 src/main.c src/${INCLUDE_LIBRARY_C} include/${INCLUDE_LIBRARY_H})
add_executable(${PROJECT_NAME}_O3 src/main.c src/${INCLUDE_LIBRARY_C} include/${INCLUDE_LIBRARY_H})
add_executable(${PROJECT_NAME}_seq_O1 src/sequential.c src/${INCLUDE_LIBRARY_SEQ_C} include/${INCLUDE_LIBRARY_SEQ_H})
add_executable(${PROJECT_NAME}_seq_O2 src/sequential.c src/${INCLUDE_LIBRARY_SEQ_C} include/${INCLUDE_LIBRARY_SEQ_H})
add_executable(${PROJECT_NAME}_seq_O3 src/sequential.c src/${INCLUDE_LIBRARY_SEQ_C} include/${INCLUDE_LIBRARY_SEQ_H})

add_executable(WriteVectOnFile src/WriteVectOnFile.c)



target_include_directories(${PROJECT_NAME}_O1 PRIVATE include)
target_include_directories(${PROJECT_NAME}_O2 PRIVATE include)
target_include_directories(${PROJECT_NAME}_O3 PRIVATE include)
target_include_directories(${PROJECT_NAME}_seq_O1 PRIVATE include)
target_include_directories(${PROJECT_NAME}_seq_O2 PRIVATE include)
target_include_directories(${PROJECT_NAME}_seq_O3 PRIVATE include)

target_include_directories(WriteVectOnFile PRIVATE include)

#link OpenMP library for parallel program

find_package(MPI REQUIRED)
if (MPI_C_FOUND)
    message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")
	target_link_libraries(${PROJECT_NAME}_O1 PUBLIC MPI::MPI_C)
    target_link_libraries(${PROJECT_NAME}_O2 PUBLIC MPI::MPI_C)
    target_link_libraries(${PROJECT_NAME}_O3 PUBLIC MPI::MPI_C)
endif()


#Set compile opt option.
target_compile_options(${PROJECT_NAME}_O1 PRIVATE -O1)
target_compile_options(${PROJECT_NAME}_O2 PRIVATE -O2)
target_compile_options(${PROJECT_NAME}_O3 PRIVATE -O3)
target_compile_options(${PROJECT_NAME}_seq_O1 PRIVATE -O1)
target_compile_options(${PROJECT_NAME}_seq_O2 PRIVATE -O2)
target_compile_options(${PROJECT_NAME}_seq_O3 PRIVATE -O3)



#set compiler features required when compiling a given target!
#Ci dice cosa ci serve per compilare codesto target
target_compile_features(${PROJECT_NAME}_O1 PRIVATE c_std_11)
target_compile_features(${PROJECT_NAME}_O2 PRIVATE c_std_11)
target_compile_features(${PROJECT_NAME}_O3 PRIVATE c_std_11)
target_compile_features(${PROJECT_NAME}_seq_O1 PRIVATE c_std_11)
target_compile_features(${PROJECT_NAME}_seq_O2 PRIVATE c_std_11)
target_compile_features(${PROJECT_NAME}_seq_O3 PRIVATE c_std_11)
target_compile_features(WriteVectOnFile PRIVATE c_std_11)


# ---------------------------------------- TEST -----------------------------------------
add_executable(my_test test/testRadixSort.c src/${INCLUDE_LIBRARY_C} include/${INCLUDE_LIBRARY_H})
target_include_directories(my_test PRIVATE include)
target_link_libraries(my_test PUBLIC MPI::MPI_C)
set (my_test mpirun -np 2 "./${name}")

add_executable(my_test2 test/testReadWriteFile.c src/${INCLUDE_LIBRARY_C} include/${INCLUDE_LIBRARY_H})
target_include_directories(my_test2 PRIVATE include)
target_link_libraries(my_test2 PUBLIC MPI::MPI_C)
set (my_test2 mpirun -np 2 "./${name}")

add_test (NAME my_test COMMAND my_test )
add_test (NAME my_test2 COMMAND my_test2  )



# ------------------------------------ DOCUMENTATION ------------------------------------
#It's generates docs for our project
find_package(Doxygen
            REQUIRED dot)
#All file needs the docs 
if(DOXYGEN_FOUND)
	set(DOXYGEN_GENERATE_XML YES)
	doxygen_add_docs(
		docs
		./src/${INCLUDE_LIBRARY_C}
		./src/${INCLUDE_LIBRARY_SEQ_C}
		)
endif()

# ---------------------------------- COMMAND  ----------------------------------
#With this command we execute the bash file
add_custom_target (
    first_command
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/measurements/analysis.bash ${CMAKE_CURRENT_BINARY_DIR} ${PROJECT_NAME}
    COMMENT "FIRST COMMAND EXECUTION [${CMAKE_CURRENT_SOURCE_DIR}/measurements/analisis.bash ${CMAKE_CURRENT_BINARY_DIR}] ${PROJECT_NAME}"
    DEPENDS	
        ${PROJECT_NAME}_O1
        ${PROJECT_NAME}_O2
        ${PROJECT_NAME}_O3
        ${PROJECT_NAME}_seq_O1
        ${PROJECT_NAME}_seq_O2
        ${PROJECT_NAME}_seq_O3
        WriteVectOnFile
)

# --------------------------------- OTHER COMMAND  ---------------------------------
#Python3 for extract the data from analisis, PipEnv because we like it!
find_package(Python3 COMPONENTS Interpreter REQUIRED)
find_program(PipEnv pipenv REQUIRED)

add_custom_target(
	install_pipenv
	COMMAND PIPENV_PIPFILE=${CMAKE_CURRENT_SOURCE_DIR}/Pipfile ${PipEnv} install --ignore-pipfile
    COMMENT "Installing pipenv dependencies [${PipEnv} install]"
    DEPENDS	${CMAKE_CURRENT_SOURCE_DIR}/Pipfile.lock
)

add_custom_target(
    second_command
    COMMAND PIPENV_PIPFILE=${CMAKE_CURRENT_SOURCE_DIR}/Pipfile ${PipEnv} run python ${CMAKE_CURRENT_SOURCE_DIR}/measurements/extract.py
    COMMENT "Generating version header [${PipEnv} run python ${CMAKE_CURRENT_SOURCE_DIR}/measures/extract.py]"
    DEPENDS	
		install_pipenv
		${CMAKE_CURRENT_SOURCE_DIR}/measurements/extract.py
)


