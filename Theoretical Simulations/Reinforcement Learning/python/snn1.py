import snntorch as snn
import torch
from snntorch import utils
from torchvision import datasets, transforms
from torch.utils.data import DataLoader

batch_size = 128
data_path = '/tmp/data/mnist'
num_classes = 10

dtype = torch.float

transform = transforms.Compose([
    transforms.Resize((28,28)),
    transforms.Grayscale(),
    transforms.ToTensor(),
    transforms.Normalize((0,),(1,))
])

mnist_train = datasets.MNIST(data_path,train=True,download=True,transform=transform)

subset = 10
mnist_train = utils.data_subset(mnist_train,subset)

train_loader = DataLoader(mnist_train,batch_size=batch_size,shuffle=True)
#   Spike Encoding
#   Rate Coding of MNIST
num_steps = 100
raw_vector = torch.ones(num_steps)*0.5
rate_coded_vector = torch.bernoulli(raw_vector)

print(f"The output is spiking {rate_coded_vector.sum()*100/len(rate_coded_vector):.2f}% of the time.")