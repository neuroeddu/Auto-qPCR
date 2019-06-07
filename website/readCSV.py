import csv

def csvreader():
    with open('test.csv') as csv_file:
        csv_reader = csv.reader(csv_file)
        for row in csv_reader:
            print(",".join(row))


def csvwriter(stream, str):
    stream = stream + "\n " + str
    return stream.replace("=", ",")


def csvwriter2(stream):
    stream = stream + "\n Goodbye, world"
    return stream.replace("=", ",")