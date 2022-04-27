# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of sideloaders
"""

import os
from typing import Any, Dict, List, Optional

from antismash.common.module_results import DetectionResults
from antismash.common.secmet.record import Record
from antismash.config import ConfigType, get_config
from antismash.config.args import ModuleArgs

from .mibig import mibig_loader, MibigAnnotations
from .html_output import will_handle, generate_html

NAME = "mibig"
SHORT_DESCRIPTION = "module for MIBiG mode"

def get_arguments() -> ModuleArgs:
    """ Constructs commandline arguments and options for this module
    """
    args = ModuleArgs("MIBiG Mode", "mibig")
    args.add_option('--mibig-mode',
                    dest='mibig_mode',
                    action='store_true',
                    default=False,
                    help="Generate MiBIG output instead of running all antiSMASH analyses")
    args.add_option("mibig-json",
                    dest="mibig_json",
                    metavar="JSON",
                    type=str,
                    default="",
                    help=("Sideload MIBiG annotations from the JSON file in the given path."))
    return args


def check_options(options: ConfigType) -> List[str]:
    """ Checks the options to see if there are any issues before
        running any analyses
    """
    if options.mibig_json:
        if not os.path.exists(options.mibig_json):
            return ["MIBiG annotation JSON cannot be found at '%s'" % options.mibig_json]
    return []


def is_enabled(options: ConfigType) -> bool:
    """  Uses the supplied options to determine if the module should be run
    """
    return options.mibig_mode and options.mibig_json


def regenerate_previous_results(results: Dict[str, Any], record: Record,
                                options: ConfigType) -> Optional[MibigAnnotations]:
    """ Regenerate previous results. """
    return None # unimplemented


def run_on_record(record: Record, previous_results: Optional[MibigAnnotations],
                  options: ConfigType) -> DetectionResults:
    """ Annotate record using MIBiG JSON.
    """
    if previous_results:
        return previous_results
    return mibig_loader(options.mibig_json, record)


def prepare_data(_logging_only: bool = False) -> List[str]:
    """ Ensures packaged data is fully prepared

        Arguments:
            logging_only: whether to return error messages instead of raising exceptions

        Returns:
            a list of error messages (only if logging_only is True)
    """
    return []


def check_prereqs(options: ConfigType) -> List[str]:
    """ Check that prerequisites are satisfied.
    """
    return []
