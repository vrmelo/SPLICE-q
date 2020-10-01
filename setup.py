import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(name='SPLICE-q',
                 version='1.0.0',
                 scripts=["spliceq/SPLICE-q.py"],
                 description='SPLICE-q is a fast and user-friendly Python tool for genome-wide SPLICing Efficiency quantification',
                 long_description=long_description,
                 long_description_content_type="text/markdown",
                 url='https://github.com/vrmelo/SPLICE-q',
                 author='Verônica R Melo Costa, Julianus Pfeuffer, Annita Louloupi, Ulf A V Ørom, Rosario M Piro',
                 author_email='veronica.melocosta@gmail.com',
                 license='GPL-2',
                 packages=setuptools.find_packages(),
                 zip_safe=False,
                 classifiers=[
                     "Programming Language :: Python :: 3",
                     "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
                     "Operating System :: OS Independent",
                 ],
                 python_requires='>=3.6',
                 # Project uses reStructuredText, so ensure that the docutils get
                 # installed or upgraded on the target machine
                 install_requires=["pysam", "rich", "interlap", "numpy"]
                 )
