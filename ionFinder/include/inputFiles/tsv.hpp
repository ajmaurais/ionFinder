//
// tsv.hpp
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

#ifndef tsv_hpp
#define tsv_hpp

#include <inputFiles/inputFiles.hpp>
#include <scanData.hpp>
#include <utils.hpp>
#include <tsvFile.hpp>
#include <tsv_constants.hpp>

namespace inputFiles {
	const std::string TSV_INPUT_REQUIRED_COLNAMES [] = {SAMPLE_NAME, SEQUENCE, PRECURSOR_FILE, SCAN_NUM};
	int const TSV_INPUT_REQUIRED_COLNAMES_LEN = 4;
	const std::string TSV_INPUT_OPTIONAL_COLNAMES [] = {PARENT_ID, PARENT_PROTEIN, PARENT_DESCRIPTION, MATCH_DIRECTION,
                                                        FORMULA, FULL_SEQUENCE, UNIQUE, CHARGE, SCORE, PRECURSOR_MZ,
                                                        PRECURSOR_SCAN};
	int const TSV_INPUT_OPTIONAL_COLNAMES_LEN = 10;

	class Tsv : public InputFile {
	private:
        bool read(const std::string &ifname, std::vector<inputFiles::Scan> &scans,
                  const std::string& sampleName) const override;
    public:
	    Tsv() : InputFile() {
	        _fileType = InputFileType::TSV;
	        _fileExtension = "tsv";
	    }

        bool findInputFiles(const std::vector<std::string>& inputArgs, std::string& wd) override;
    };
}

#endif /* tsv_hpp */
