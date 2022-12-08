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
#include <Evaluator/detail/Evaluator.h>

namespace jana::parser {
  
  // parsing result template
  template<class ResultType>
    struct Result {
      bool success;
      ResultType result;
    };

  class Parser {
    public:

      Parser(bool debug=false, int min_precision=6) : m_debug(debug), m_min_precision(min_precision) {

        //
        // FIXME: decide which unit system to use
        //

        // Geant4 unit system 
        // 1 = millimeter = MeV = nanosecond = positron charge = kelvin = mole = candela = radian = steradian
        /*
        m_eval = std::make_unique<dd4hep::tools::Evaluator::Object>(
            1.e+3, 1./1.60217733e-25, 1.e+9, 1./1.60217733e-10, 1.0, 1.0, 1.0);
            */

        // Geant4 convention, using DD4hep names
        /*
        m_eval = std::make_unique<dd4hep::tools::Evaluator::Object>(
            dd4hep::millimeter,
            dd4hep::gram * 1e-6,
            dd4hep::nanosecond,
            dd4hep::ampere,
            dd4hep::kelvin,
            dd4hep::mole,
            dd4hep::candela,
            dd4hep::radian
            );
            */

        // DD4hep default system
        // 1 = centimeter = GeV = second = nanoampere = kelvin = mole = candela = radian
        // /*
        m_eval = std::make_unique<dd4hep::tools::Evaluator::Object>(
            dd4hep::meter,
            dd4hep::kilogram,
            dd4hep::second,
            dd4hep::ampere,
            dd4hep::kelvin,
            dd4hep::mole,
            dd4hep::candela,
            dd4hep::radian
            );
            // */
      };
      ~Parser() {};

      // evaluate DD4hep expression `expr` and return a string
      // - if `expr` is itself a string, return `expr` silently
      // - this is useful for parsing CLI options
      // - optionally specify `key` to clarify debugging printouts
      Result<std::string> dd4hep_to_string(const std::string& expr, const std::string& key="") {
        if(m_debug) fmt::print("PARSE: dd4hep_to_string({}){}\n", expr, key==""? "" : " for key "+key);

        // handle empty string, which would otherwise cause en error in parsing
        if(expr=="") return { true, expr };

        // call Evaluator::evaluate to parse units and do the math
        auto parsed = m_eval->evaluate(expr.c_str());

        // return a `Result`, with the calculated number re-stringified without any units
        if(m_debug) fmt::print("-> status: {}\n", parsed.status());
        switch(parsed.status()) {

          case dd4hep::tools::Evaluator::OK:
            { // likely a number that was parsed successfuly; stringify it
              std::stringstream ss;
              // increase stringify precision, if the user specified higher than `m_min_precision`
              auto precision = std::max(
                  (unsigned long) m_min_precision,
                  expr.substr(0,expr.find("*")).size() // count `char`s up to first `*`
                  );
              // stringify
              ss << std::setprecision(precision) << parsed.result();
              auto result = ss.str();
              if(m_debug) fmt::print("-> evaluator result: {}\n-> string: '{}'\n", parsed.result(), result);
              return { true, result };
            }

          case dd4hep::tools::Evaluator::ERROR_UNKNOWN_VARIABLE:
            { // likely a string (as long as it's not any specific unit name); return `expr` as is
              if(expr.find('*') != std::string::npos) // try to detect units typos
                fmt::print(stderr,"WARNING: parsing '{}' as a string; is there a typo in the units?\n",expr);
              if(m_debug) fmt::print("-> string: '{}'\n", expr);
              return { true, expr };
            }

          case dd4hep::tools::Evaluator::ERROR_UNEXPECTED_SYMBOL:
            { // unexpected symbol; if it's just a comma, attempt to parse as a list, otherwise complain and return `expr` as is
              if(expr.find(',') != std::string::npos) {
                std::string result = "";
                bool success       = true;
                if(m_debug) fmt::print("{:-^30}\n"," begin list ");
                // tokenize
                std::istringstream expr_s(expr);
                std::string tok;
                while(getline(expr_s, tok, ',')) { // call `dd4hep_to_string` on each list element
                  auto parsed_tok = dd4hep_to_string(tok);
                  result += "," + parsed_tok.result;
                  success &= parsed_tok.success; // innocent until proven guilty
                }
                result.erase(0,1); // remove leading comma
                if(m_debug) fmt::print("{:-^30}\n-> string: {}\n", " end list ", result);
                return { success, result };
              }
              break;
            }
        };

        // likely an error; complain and return `expr` as is
        fmt::print(stderr,"ERROR: cannot evaluate '{}'; ",expr);
        parsed.print_error();
        return { false, expr };
      }

    private:
      std::unique_ptr<dd4hep::tools::Evaluator::Object> m_eval;
      bool m_debug;
      int m_min_precision;
  };
}
