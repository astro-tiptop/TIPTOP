import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, TensorDataset, random_split, DataLoader 

#from torch_geometric.nn import GCNConv, global_mean_pool, global_add_pool
#from torch_geometric.data import Data
#from torch_geometric.loader import DataLoader

import torch
import torch.nn as nn
from torchmetrics.regression import SymmetricMeanAbsolutePercentageError
from torchmetrics.regression import MeanAbsolutePercentageError
from torchmetrics.regression import MeanAbsoluteError


device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
#device = torch.device("cpu")
# torch.set_default_dtype(torch.float64)

class TriangleDataset(Dataset):
    def __init__(self, inputs, targets):
        self.inputs = inputs
        self.targets = targets

    def __len__(self):
        return len(self.inputs)

    def __getitem__(self, idx):
        sample = {'input': self.inputs[idx], 'target': self.targets[idx]}
        return sample

def getLoader(inputs, targets, batch_size):    
    dataset = TriangleDataset(inputs, targets)
    # DataLoader
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=True, drop_last=True)
    return loader

class NeuralNetwork(nn.Module):
    def __init__(self, ninputs, l_sizes):
        super(NeuralNetwork, self).__init__()
        self.nl = len(l_sizes)

        self.layersA1 = nn.ModuleList()
        self.layersA1.append(nn.Linear(4, l_sizes[0]))
        for s in range(1, self.nl):
            self.layersA1.append(nn.Linear(l_sizes[s-1], l_sizes[s]))
        self.layersA2 = nn.ModuleList()
        self.layersA2.append(nn.Linear(4, l_sizes[0]))
        for s in range(1, self.nl):
            self.layersA2.append(nn.Linear(l_sizes[s-1], l_sizes[s]))
        self.layersA3 = nn.ModuleList()
        self.layersA3.append(nn.Linear(4, l_sizes[0]))
        for s in range(1, self.nl):
            self.layersA3.append(nn.Linear(l_sizes[s-1], l_sizes[s]))

        self.layersB1 = nn.ModuleList()
        self.layersB1.append(nn.Linear(3, l_sizes[0]))
        for s in range(1, self.nl):
            self.layersB1.append(nn.Linear(l_sizes[s-1], l_sizes[s]))
        self.layersB2 = nn.ModuleList()
        self.layersB2.append(nn.Linear(3, l_sizes[0]))
        for s in range(1, self.nl):
            self.layersB2.append(nn.Linear(l_sizes[s-1], l_sizes[s]))
        self.layersB3 = nn.ModuleList()
        self.layersB3.append(nn.Linear(3, l_sizes[0]))
        for s in range(1, self.nl):
            self.layersB3.append(nn.Linear(l_sizes[s-1], l_sizes[s]))
        self.layersB4 = nn.ModuleList()
        self.layersB4.append(nn.Linear(3, l_sizes[0]))
        for s in range(1, self.nl):
            self.layersB4.append(nn.Linear(l_sizes[s-1], l_sizes[s]))

        self.layersC = nn.ModuleList()
        self.layersC.append(nn.Linear(7*l_sizes[0], 7*l_sizes[1]))
        for s in range(1, self.nl-1):
            self.layersC.append(nn.Linear(7*l_sizes[s], 7*l_sizes[s+1]))
        self.layersC.append(nn.Linear(7*l_sizes[-1], 1))
                
        self.dp = nn.Dropout(p=0.01)
        self.leaky_relu1 = nn.LeakyReLU(0.02)
        self.sigma = nn.Sigmoid()
        self.elu = nn.ELU(0.02)
        self.tan = nn.Tanh()
        self.act1 = self.elu # self.leaky_relu1 # self.leaky_relu1
        self.act2 = self.elu # self.tan        
        self.act3 = self.elu
        self.double()
        
    def forward(self, xy):
        x1 = xy[:, 0::3]
        x1 = self.layersA1[0](x1)
        for l in range(1,self.nl):
            x1 = self.dp(self.act1(self.layersA1[l](x1)))
        x2 = xy[:, 1::3]
        x2 = self.layersA2[0](x2)
        for l in range(1,self.nl):
            x2 = self.dp(self.act1(self.layersA2[l](x2)))
        x3 = xy[:, 2::3]
        x3 = self.layersA3[0](x3)
        for l in range(1,self.nl):
            x3 = self.dp(self.act1(self.layersA3[l](x3)))       
            
        xx1 = xy[:, :3]
        xx1 = self.layersB1[0](xx1)
        for l in range(1,self.nl):
            xx1 = self.dp(self.act1(self.layersB1[l](xx1)))
        xx2 = xy[:, 3:6]
        xx2 = self.layersB2[0](xx2)
        for l in range(1,self.nl):
            xx2 = self.dp(self.act1(self.layersB2[l](xx2)))
        xx3 = xy[:, 6:9]
        xx3 = self.layersB3[0](xx3)
        for l in range(1,self.nl):
            xx3 = self.dp(self.act1(self.layersB3[l](xx3)))
        xx4 = xy[:, 9:]
        xx4 = self.layersB4[0](xx4)
        for l in range(1,self.nl):
            xx4 = self.dp(self.act1(self.layersB4[l](xx4)))

        z1 = torch.cat((xx1, xx2, xx3, xx4), 1)
#       y = xy[:, 6:]
#       y = self.layersB[0](y)
#       for l in range(1,self.nl):
#           y = self.act2(self.layersB[l](y))
#           y = self.dp(y)
#        z = torch.cat((x, y), 1)
        z2 = torch.cat((x1, x2, x3), 1)

        zz = torch.cat((z1, z2), 1)

        for l in range(self.nl-1):
            zz = self.act3(self.layersC[l](zz))
            zz = self.dp(zz)
        zz = self.layersC[-1](zz)
        return zz
        
        
    def save_model(self, path):
        """ Save the model parameters to the specified file. """
        torch.save(self.state_dict(), path)
        print(f"Model saved to {path}")

    def load_model(self, path):
        """ Load model parameters from the specified file. """
        self.load_state_dict(torch.load(path))
        self.double()
        self.eval()  # Set the model to evaluation mode
        print(f"Model loaded from {path}")

    def setData(self, inputDataT, jitterT, testShare, updateScales = True):
        with torch.no_grad():
            if updateScales:
                self.maxsJ = torch.max(jitterT, 0).values
                self.minsJ = torch.min(jitterT, 0).values
                self.maxsI = torch.max(inputDataT, 0).values
                self.minsI = torch.min(inputDataT, 0).values                
                self.maxsI[7] = self.maxsI[8] = self.maxsI[6]
                self.minsI[7] = self.minsI[8] = self.minsI[6]
#            print("jitterT self.maxs", self.maxsJ)
#            print("jitterT self.mins", self.minsJ)
#            print("inputDataT self.maxs", self.maxsI)
#            print("inputDataT self.mins", self.minsI)            
            jitterT = 2.0 * (jitterT - self.minsJ) / (self.maxsJ - self.minsJ) + 0.1
            inputDataT = (inputDataT - self.minsI) / (self.maxsI - self.minsI) - 0.1
#            inputDataT = inputDataT / self.maxsI
            
        self.inputDataT = inputDataT.to(device)
        self.jitterT = jitterT.to(device)
        self.dataset = TensorDataset(self.inputDataT, self.jitterT)        
        self.test_split = testShare
        self.batch_size = 64
        self.total_samples = len(self.dataset)
        self.test_size = int(self.test_split * self.total_samples)
        self.train_size = self.total_samples - self.test_size
        torch.manual_seed(5125)
        if testShare<1.0:
            self.train_dataset, self.test_dataset = random_split(self.dataset, [self.train_size, self.test_size])
        else:
            self.train_dataset, self.test_dataset = None, self.dataset
#        print('train_dataset', self.train_size)
#        print('test_dataset', self.test_size)
        # Creating data loaders for training and testing
        
#        if testShare<1.0:
#            self.train_loader = getLoader(self.train_dataset[:][0], self.train_dataset[:][1], self.batch_size)
#        self.test_loader = getLoader(self.test_dataset[:][0], self.test_dataset[:][1], self.batch_size)
        
# GNexp        
        if testShare<1.0:
            self.train_loader = DataLoader(self.train_dataset, batch_size=self.batch_size, shuffle=True, drop_last=True)
        self.test_loader = DataLoader(self.test_dataset, batch_size=self.batch_size, shuffle=True, drop_last=True)

    def trainModel(self, num_epochs, steps, modelNameFile):
        self.to(device)        
        bestLossTest = 1e9
        for ss in steps:
            self.train()        
            self.criterion0 = nn.MSELoss().to(device)
            self.criterion1 = MeanAbsolutePercentageError().to(device) #
            optimizer = optim.Adam(self.parameters(), lr=ss, weight_decay=0 )
            for epoch in range(num_epochs):
                
#                for batch_data in self.train_loader:
#                    inputs = batch_data['input']
#                    targets = batch_data['target']
#                    batch_size = inputs.size(0)
#                    num_nodes_per_graph = 3  # Since we have 3 vertices per triangle
#                    num_nodes = batch_size * num_nodes_per_graph
#                    # Create edge_index and batch tensors
#                    edge_index = torch.tensor([[0, 1, 2], [1, 2, 0]], dtype=torch.long).repeat(batch_size, 1, 1).view(2, -1)
#                    batch_tensor = torch.arange(batch_size).repeat_interleave(num_nodes_per_graph)
#                    optimizer.zero_grad()
#                    output = self(inputs, edge_index, batch_tensor)
#                    loss = self.criterion1(output, targets)
#                    loss.backward()
#                    optimizer.step()                    
                    
                for inputs, targets in self.train_loader:
                    self.train()
                    optimizer.zero_grad()  # Clear existing gradients
                    outputs = self(inputs)  # Forward pass
                    if epoch<5:
                        loss = self.criterion0(outputs, targets)  # Compute loss
                    else:
                        loss = self.criterion1(outputs, targets)  # Compute loss                        
                    loss.backward()  # Backward pass
                    optimizer.step()  # Update model parameters
                    
                lossTest = self.testModel()
                print(f'Loss at Epoch [{epoch+1}/{num_epochs}], Loss: {loss.item():.5f}, Loss Val: {lossTest:.5f}')
                if lossTest < bestLossTest and epoch>10:
                    bestLossTest = lossTest
                    self.save_model(modelNameFile)
                    print(f'Best Loss at Epoch [{epoch+1}/{num_epochs}], Loss: {loss.item():.5f}, Loss Val: {lossTest:.5f}')
                
        inputDataTest, jitterTest = self.test_loader.dataset[:]
        inputDataTrain, jitterTrain = self.train_loader.dataset[:]
        with torch.no_grad():
            jitterApproxTrain = self(inputDataTrain)
            jitterApproxTest = self(inputDataTest)
            
        jitterTrain = jitterTrain.detach().cpu().numpy()
        jitterApproxTrain = jitterApproxTrain.detach().cpu().numpy()
        jitterApproxTrain = jitterApproxTrain

        jitterTest = jitterTest.detach().cpu().numpy()
        jitterApproxTest = jitterApproxTest.detach().cpu().numpy()
        jitterApproxTest = jitterApproxTest
        
        return inputDataTrain.detach().cpu().numpy(), jitterTrain, jitterApproxTrain, inputDataTest.detach().cpu().numpy(), jitterTest, jitterApproxTest


    def testModel(self):
        self.eval()
        running_loss2 = 0.0
        with torch.no_grad():  # No need to track gradients during evaluation
            for inputsT, labelsT in self.test_loader:
                outputsT = self(inputsT)
                lossV = self.criterion1(outputsT, labelsT)
                running_loss2 += lossV.cpu().item()
#            for batch_data in self.test_loader:
#                inputs = batch_data['input']
#                targets = batch_data['target']
#                batch_size = inputs.size(0)
#                num_nodes_per_graph = 3  # Since we have 3 vertices per triangle
#                num_nodes = batch_size * num_nodes_per_graph
#                # Create edge_index and batch tensors
#                edge_index = torch.tensor([[0, 1, 2], [1, 2, 0]], dtype=torch.long).repeat(batch_size, 1, 1).view(2, -1)
#                batch_tensor = torch.arange(batch_size).repeat_interleave(num_nodes_per_graph)
#                output = self(inputs, edge_index, batch_tensor)
#                lossV = self.criterion1(output, targets)
#                running_loss2 += lossV.cpu().item()

        self.train()
#        print('error scale',  float((self.test_size / self.batch_size)))
        running_loss2 = running_loss2 / int(self.test_size / self.batch_size)
        return running_loss2 # lossV.item()
'''     
    def __init__(self):
        super(NeuralNetwork, self).__init__()
        self.conv1 = GCNConv(4, 256)  # Adjust input features from 2 to 4
        self.conv2 = GCNConv(256, 256)
        self.conv3 = GCNConv(256, 128)
        self.conv4 = GCNConv(128, 64)
        self.conv5 = GCNConv(64, 64)
        self.fc1 = nn.Linear(64, 64)
        self.fc2 = nn.Linear(64, 64)
        self.fc3 = nn.Linear(64, 1)
        self.relu = nn.LeakyReLU(0.1) # nn.ReLU()
        self.tan = nn.Tanh()
        self.double()

        
    def forward(self, inputs, edge_index, batch):
        # Reshape the input tensor
        x_coords = inputs[:, :3].reshape(-1, 1)
        y_coords = inputs[:, 3:6].reshape(-1, 1)
        scalar1 = inputs[:, 6:9].reshape(-1, 1)
        scalar2 = inputs[:, 9:12].reshape(-1, 1)
        x = torch.cat([x_coords, y_coords, scalar1, scalar2], dim=1)  # Shape: (batch_size*3, 4)
        x = self.conv1(x, edge_index)
        x = self.relu(self.conv2(x, edge_index))
        x = self.relu(self.conv3(x, edge_index))
        x = self.relu(self.conv4(x, edge_index))
        x = self.relu(self.conv5(x, edge_index))
        x = global_mean_pool(x, batch)  # Global pooling
        x = self.relu(self.fc1(x))
        x = self.fc2(x)
        x = self.fc3(x)
        return x

        
'''
