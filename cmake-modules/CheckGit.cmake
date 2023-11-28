# This CheckGit functionality comes in part from https://gitlab.com/jhamberg/cmake-examples
# and has been modified by David Van Komen for the purposes of this project

# ORIGINAL LICENSE:
#
# MIT License
#
# Copyright (c) 2901 Jonathan Hamberg
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#

# Allow if statements to have bare includes...
# set(CMAKE_POLICY_DEFAULT_CMP0012 NEW)

set(CURRENT_LIST_DIR ${CMAKE_CURRENT_LIST_DIR})
if (NOT DEFINED pre_configure_dir)
    set(pre_configure_dir ${CMAKE_CURRENT_LIST_DIR})
endif ()

if (NOT DEFINED post_configure_dir)
    set(post_configure_dir ${CMAKE_BINARY_DIR}/generated)
endif ()

set(pre_configure_file ${pre_configure_dir}/git_version_and_date.cpp.in)
set(post_configure_file ${post_configure_dir}/git_version_and_date.cpp)

function(CheckGitWrite git_hash)
    file(WRITE ${CMAKE_BINARY_DIR}/git-state.txt ${git_hash})
endfunction()

function(CheckGitRead git_hash)
    if (EXISTS ${CMAKE_BINARY_DIR}/git-state.txt)
        file(STRINGS ${CMAKE_BINARY_DIR}/git-state.txt CONTENT)
        LIST(GET CONTENT 0 var)

        set(${git_hash} ${var} PARENT_SCOPE)
    endif ()
endfunction()

function (CheckGitDirtyWrite dirty_status)
    file(WRITE ${CMAKE_BINARY_DIR}/dirty-status.txt ${dirty_status})
endfunction()

function(CheckGitDirtyRead dirty_status)
    if (EXISTS ${CMAKE_BINARY_DIR}/dirty-status.txt)

        file(STRINGS ${CMAKE_BINARY_DIR}/dirty-status.txt CONTENT)
        LIST(GET CONTENT 0 var)
        set(${dirty_status} ${var} PARENT_SCOPE)
    endif ()
endfunction()

function(CheckGitVersion)
    # Get the latest abbreviated commit hash of the working branch
    execute_process(
        COMMAND git log -1 --format=%h
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE
        )

    execute_process(
        COMMAND "git" "diff" "--quiet"
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        RESULT_VARIABLE GIT_DIRTY_STAT_1
        )
    execute_process(
        COMMAND "git" "diff" "--cached" "--quiet"
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        RESULT_VARIABLE GIT_DIRTY_STAT_2
        )

    if (GIT_DIRTY_STAT_1 EQUAL "1" OR GIT_DIRTY_STAT_2 EQUAL "1")
        set(GIT_DIRTY "-dirty")
    else()
        set(GIT_DIRTY " (clean repo)")
    endif()

    execute_process(
        COMMAND date
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE COMPILE_DATE_TIME
        OUTPUT_STRIP_TRAILING_WHITESPACE
        )

    ### ---- CHECK IF GIT HASH IS CACHED
    CheckGitRead(GIT_HASH_CACHE)
    if (NOT EXISTS ${post_configure_dir})
        file(MAKE_DIRECTORY ${post_configure_dir})
    endif ()

    if (NOT EXISTS ${post_configure_dir}/git_version_and_date.h)
        file(COPY ${pre_configure_dir}/git_version_and_date.h DESTINATION ${post_configure_dir})
    endif()

    if (NOT DEFINED GIT_HASH_CACHE)
        set(GIT_HASH_CACHE "INVALID")
    endif ()

    set(please_configure OFF)

    # Only update the git_version.cpp if the hash has changed. This will
    # prevent us from rebuilding the project more than we need to.
    if (NOT ${GIT_HASH} STREQUAL ${GIT_HASH_CACHE} OR NOT EXISTS ${post_configure_file})
        # Set the GIT_HASH_CACHE variable the next build won't have
        # to regenerate the source file.
        CheckGitWrite("${GIT_HASH}")

        set(please_configure ON)

    endif ()
    ### ---- END GIT HASH

    ### ---- CHECK GIT DIRTY STATUS
    CheckGitDirtyRead(GIT_DIRTY_CACHE)

    if (NOT DEFINED GIT_DIRTY_CACHE)
        set(GIT_DIRTY_CACHE "INVALID")
    endif ()

    if (NOT ${GIT_DIRTY} STREQUAL ${GIT_DIRTY_CACHE} OR NOT EXISTS ${post_configure_file})

        CheckGitDirtyWrite("${GIT_DIRTY}")

        set(please_configure ON)

    endif ()
    ### ---- END GIT DIRTY STATUS

    # only configure if something has changed
    # NOTE: even if nothing has changed in the repo, every time make is called it'll update
    if (please_configure)
        configure_file(${pre_configure_file} ${post_configure_file} @ONLY)
    endif ()


endfunction()

function(CheckGitSetup)

    add_custom_target(AlwaysCheckGit COMMAND ${CMAKE_COMMAND}
        -DRUN_CHECK_GIT_VERSION=1
        -Dpre_configure_dir=${pre_configure_dir}
        -Dpost_configure_file=${post_configure_dir}
        -DGIT_HASH_CACHE=${GIT_HASH_CACHE}
        -P ${CURRENT_LIST_DIR}/CheckGit.cmake
        BYPRODUCTS ${post_configure_file}
        )

    add_library(dendro_git_version_and_date ${CMAKE_BINARY_DIR}/generated/git_version_and_date.cpp)
    target_include_directories(dendro_git_version_and_date PUBLIC ${CMAKE_BINARY_DIR}/generated)
    add_dependencies(dendro_git_version_and_date AlwaysCheckGit)

    CheckGitVersion()
endfunction()

# This is used to run this function from an external cmake process.
if (RUN_CHECK_GIT_VERSION)
    CheckGitVersion()
endif ()

