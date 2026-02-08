import torch
import torch.nn as nn


class ResidualBlock_v2(nn.Module):
    def __init__(self, channels, width, dilation):
        super(ResidualBlock_v2, self).__init__()
        self.bn1 = nn.BatchNorm1d(channels)
        self.relu = nn.ReLU()
        self.conv1 = nn.Conv1d(
            channels, channels, width, dilation=dilation, padding="same"
        )
        self.bn2 = nn.BatchNorm1d(channels)
        self.conv2 = nn.Conv1d(
            channels, channels, width, dilation=dilation, padding="same"
        )

    def forward(self, x):
        residual = x
        out = self.bn1(x)
        out = self.relu(out)
        out = self.conv1(out)
        out = self.bn2(out)
        out = self.relu(out)
        out = self.conv2(out)
        return out + residual


class TRANSAID_Embedding(nn.Module):
    def __init__(self, embedding_dim=128, output_channels=3):
        super(TRANSAID_Embedding, self).__init__()

        # Embedding层：5个可能的值(0-4)，0用作padding
        self.embedding = nn.Embedding(
            num_embeddings=5,  # 0=padding, 1=A, 2=C, 3=G, 4=T/U
            embedding_dim=embedding_dim,
            padding_idx=0,  # 指定padding_idx确保padding token的embedding始终为0
        )

        self.conv1 = nn.Conv1d(embedding_dim, 32, kernel_size=3, padding="same")
        self.bn1 = nn.BatchNorm1d(32)
        self.relu = nn.ReLU()

        self.rb1 = nn.Sequential(*[ResidualBlock_v2(32, 26, 1) for _ in range(4)])
        self.conv2 = nn.Conv1d(32, 64, kernel_size=1, padding="same")
        self.bn2 = nn.BatchNorm1d(64)

        self.rb2 = nn.Sequential(*[ResidualBlock_v2(64, 26, 2) for _ in range(4)])
        self.conv3 = nn.Conv1d(64, 128, kernel_size=1, padding="same")
        self.bn3 = nn.BatchNorm1d(128)

        self.rb3 = nn.Sequential(*[ResidualBlock_v2(128, 36, 5) for _ in range(4)])
        self.conv4 = nn.Conv1d(128, 32, kernel_size=1, padding="same")

        self.conv5 = nn.Conv1d(32, output_channels, kernel_size=1, padding="same")

    def forward(self, x):
        # x形状: (batch_size, seq_len)，值域: 0-4
        x = self.embedding(x)  # (batch_size, seq_len, embedding_dim)
        x = x.permute(0, 2, 1)  # (batch_size, embedding_dim, seq_len)

        x = self.relu(self.bn1(self.conv1(x)))
        x = self.rb1(x)
        x = self.relu(self.bn2(self.conv2(x)))
        x = self.rb2(x)
        x = self.relu(self.bn3(self.conv3(x)))
        x = self.rb3(x)
        x = self.conv4(x)
        x = self.conv5(x)

        x = x.permute(0, 2, 1)  # (batch_size, seq_len, output_channels)
        return x
