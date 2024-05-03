import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset

import torch
import torch.nn as nn

class NeuralNetwork(nn.Module):
    def __init__(self, ninputs):
        super(NeuralNetwork, self).__init__()
        self.layer1 = nn.Linear(ninputs, 100)
        self.layer2 = nn.Linear(100, 200)
        self.layer3 = nn.Linear(200, 400)
        self.layer4 = nn.Linear(400, 200)
        self.layer5 = nn.Linear(200, 100)
        self.layer6 = nn.Linear(100, 50)
        self.output_layer = nn.Linear(50, 1)
        self.leaky_relu = nn.LeakyReLU(0.01)

    def forward(self, x):
        x = self.layer1(x)
        x = self.leaky_relu(self.layer2(x))
        x = self.leaky_relu(self.layer3(x))
        x = self.leaky_relu(self.layer4(x))
        x = self.leaky_relu(self.layer5(x))
        x = self.layer6(x)
        x = self.output_layer(x)
        return x

    def save_model(self, path):
        """ Save the model parameters to the specified file. """
        torch.save(self.state_dict(), path)
        print(f"Model saved to {path}")

    def load_model(self, path):
        """ Load model parameters from the specified file. """
        self.load_state_dict(torch.load(path))
        self.eval()  # Set the model to evaluation mode
        print(f"Model loaded from {path}")
