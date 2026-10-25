# 01.common.smk
# Import common utilities, no business rules here.

import os
import sys

sys.path.insert(0, os.path.join(workflow.basedir, "rules", "utils"))
from common import (
    get_platform,
    get_mixcr_preset,
    get_sample_path,
    get_sample_index,
    module_enabled,
    parameter_enabled,
)
from tools import (
    get_software,
    mixcr_analyze_cmd,
    mixcr_postanalysis_cmd,
    mixcr_downsample_cmd,
    mixcr_export_clones_cmd,
    vdjtools_convert_cmd,
    vdjtools_diversity_cmd,
    cutadapt_cmd,
)
