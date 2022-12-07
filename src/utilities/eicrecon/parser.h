// Copyright 2022, Christopher Dilks
// Subject to the terms in the LICENSE file found in the top-level directory.
//
//

#pragma once

#include <string>
#include <fmt/format.h>

// DD4Hep
#include <DD4hep/DD4hepUnits.h>
#include <Evaluator/Evaluator.h>

namespace jana::parser {
  
  // parsing result template
  template<class ResultType>
    struct Result {
      bool success;
      ResultType result;
    };

  class Parser {
    public:

      Parser() : debug(true) {
        // FIXME: how to set the 'correct' units
        // m_eval = std::make_unique<dd4hep::tools::Evaluator>(
        //     1.e+3, 1./1.60217733e-25, 1.e+9, 1./1.60217733e-10, 1.0, 1.0, 1.0); // Geant4
        m_eval = std::make_unique<dd4hep::tools::Evaluator>(
            dd4hep::meter,
            dd4hep::kilogram,
            dd4hep::second,
            dd4hep::ampere,
            dd4hep::kelvin,
            dd4hep::mole,
            dd4hep::candela,
            dd4hep::radian
            ); // try to use what DD4hep is using
      };
      ~Parser() {};

      // evaluate DD4hep expression `expr` and return a string
      // - if `expr` is itself a string, return `expr` silently
      // - this is useful for parsing CLI options
      Result<std::string> dd4hep_to_string(const std::string& expr) {
        if(debug) fmt::print("PARSE: dd4hep_to_string({})\n",expr);
        auto parsed = m_eval->evaluate(expr);
        switch(parsed.first) {
          case dd4hep::tools::Evaluator::OK: // likely a number that was parsed successfuly; stringify it
            if(debug) fmt::print("-> result: {} -> string: '{}'\n", parsed.second, std::to_string(parsed.second));
            return { true, std::to_string(parsed.second) };
          case dd4hep::tools::Evaluator::ERROR_UNKNOWN_VARIABLE: // likely a string; return `expr` as is
            if(expr.find('*') != std::string::npos) // try to detect units typos
              fmt::print(stderr,"WARNING: parsing '{}' as a string; is there a typo in the units?\n",expr);
            if(debug) fmt::print("-> string\n");
            return { true, expr };
          default: // likely an error; complain and return `expr` as is
            fmt::print(stderr,"ERROR: cannot evaluate '{}'\n",expr);
            return { false, expr };
        };
      }

    private:
      std::unique_ptr<dd4hep::tools::Evaluator> m_eval;
      bool debug;
  };
}
