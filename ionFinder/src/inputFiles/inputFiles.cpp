//
// inputFiles.cpp
// ionFinder
// -----------------------------------------------------------------------------
// MIT License
// Copyright 2022 Aaron Maurais
// -----------------------------------------------------------------------------
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
// -----------------------------------------------------------------------------
//

#include <inputFiles/inputFiles.hpp>

const std::string inputFiles::Scan::REVERSE_MATCH = "reverse_";

inputFiles::Scan& inputFiles::Scan::operator = (const inputFiles::Scan& rhs)
{
    scanData::Scan::operator=(rhs);
    _parentProtein = rhs._parentProtein;
    _parentID = rhs._parentID;
    _parentDescription = rhs._parentDescription;
    _matchDirection = rhs._matchDirection;
    _sampleName = rhs._sampleName;
    _unique = rhs._unique;

    return *this;
}

inputFiles::ModFilter inputFiles::intToModFilter(int i){
    switch(i) {
        case 0:
            return inputFiles::ModFilter::ONLY_MODIFIED;
        case 1:
            return inputFiles::ModFilter::ALL;
        case 2:
            return inputFiles::ModFilter::EXCLUDE_MODIFIED;
        default:
            throw std::invalid_argument("Can not convert '" + std::to_string(i) + "' to ModFilter!");
    }
}

/**
 \brief Get protein match direction. <br>

 Match direction is determined by checking  if \p str contains inputFiles::REVERSE_MATCH.
 \param str \p db tag from fasta header line.
 See https://www.uniprot.org/help/fasta-headers for detials on fasta headers.
 \return match direction for fasta entry
 */
inputFiles::Scan::MatchDirection inputFiles::Scan::strToMatchDirection(std::string str)
{
    if(utils::strContains(REVERSE_MATCH, utils::toLower(str)))
        return inputFiles::Scan::MatchDirection::REVERSE;
    return inputFiles::Scan::MatchDirection::FORWARD;
}
