cmake_minimum_required(VERSION 2.8.9)

project(Adquire_Data Data)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "http://www.example.com/Slicer/Extensions/Adquire Data")
set(EXTENSION_CATEGORY "Adquire_Data")
set(EXTENSION_CONTRIBUTORS "Davide Scorza (Vicomtech (San Sebastian)), Camilo Cortés (Vicomtech (San Sebastian))")
set(EXTENSION_DESCRIPTION "Adquire data")
set(EXTENSION_ICONURL "http://www.example.com/Slicer/Extensions/Adquire Data.png")
set(EXTENSION_SCREENSHOTURLS "http://www.example.com/Slicer/Extensions/Adquire Data/Screenshots/1.png")

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(CalculateDistances)
add_subdirectory(CompareDistances)
add_subdirectory(SaveFiducialList)
add_subdirectory(vtkModel)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_CPACK})
