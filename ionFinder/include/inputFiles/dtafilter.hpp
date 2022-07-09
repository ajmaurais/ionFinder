//
// dtafilter.hpp
// ionFinder
// -----------------------------------------------------------------------------
// MIT License
// Copyright 2020 Aaron Maurais
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

#ifndef dtafilter_hpp
#define dtafilter_hpp

#include <iostream>
#include <fstream>
#include <map>

#include <paramsBase.hpp>
#include <inputFiles/inputFiles.hpp>
#include <utils.hpp>

namespace inputFiles {

    bool readFilterFiles(const std::map<std::string, std::string>& filterFiles,
                         std::vector<inputFiles::Scan>&,
                         bool skipReverse,
                         ModFilter modFilter);

    class DtaFilterFile {
	private:
        static bool parse_matchDir_ID_Protein(const std::string &str, Scan &scan) ;
	public:
	    DtaFilterFile() = default;

        bool readFilterFile(const std::string& fname, const std::string& sampleName,
                                   std::vector<inputFiles::Scan>& scans,
                                   bool skipReverse = false, ModFilter modFilter = ModFilter::ALL) const;

        static void initilizeFromLine(std::string line, Scan &scan);
    };
}

#endif /* dtafilter_hpp */
