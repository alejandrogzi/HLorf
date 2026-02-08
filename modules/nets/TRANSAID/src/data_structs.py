# src/data_structs.py

from dataclasses import dataclass
from typing import List, Dict, Optional
import numpy as np


@dataclass
class SequenceFeatures:
    """序列特征数据类"""

    sequence_id: str
    sequence: str
    length: int
    gc_content: float
    features: Optional[Dict] = None


@dataclass
class ModelPrediction:
    """模型预测结果数据类"""

    sequence_id: str
    tis_probs: np.ndarray
    tts_probs: np.ndarray
    kozak_probs: np.ndarray


@dataclass
class ValidationMetrics:
    """验证指标数据类"""

    accuracy: float
    precision: float
    recall: float
    f1_score: float
    auc_score: float


@dataclass
class TranslationResult:
    """
    用于存储单个预测出的开放阅读框 (ORF) 的所有信息的标准数据结构。
    """

    # 基本信息
    sequence_id: str
    start_position: int  # 0-based index of the first base of the start codon (TIS)
    stop_position: int  # 0-based index of the first base of the stop codon (TTS)

    # 序列信息
    protein_sequence: str

    # 模型和生物学特征得分
    tis_score: float = 0.0
    tts_score: float = 0.0
    kozak_score: float = 0.0
    cai_score: float = 0.0
    gc_score: float = 0.0
    integrated_score: float = 0.0

    # 过滤状态
    passed_filter: bool = True  # Default to True, can be set to False later
    filter_reason: str = ""

    # 为了让dataclass可比较 (例如在'best'模式下), 我们让它可以被哈希
    def __hash__(self):
        return hash((self.sequence_id, self.start, self.stop))

    def __eq__(self, other):
        if not isinstance(other, TranslationResult):
            return NotImplemented
        return (self.sequence_id, self.start, self.stop) == (
            other.sequence_id,
            other.start,
            other.stop,
        )
