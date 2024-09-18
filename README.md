# JJCorrFitter
Correlation function fitting library


## Idea 1 - Inheritance Based
```mermaid

classDiagram
    Fitter1D <|-- FitterImpl
    Fitter3D <|-- FitterImpl
    CorrelationFunctionImpl <-- SourceFunctionImpl
    CorrelationFunction1D <-- InteractionTermImpl
    CorrelationFunction3D <-- InteractionTermImpl
    InteractionTermSchrodinger <|-- InteractionTermImpl
    CorrelationFunction1D <-- IntegratorOneDim
    CorrelationFunction3D <-- IntegratorMultiDim
    CorrelationFunction1D <|-- CorrelationFunctionImpl
    CorrelationFunction3D <|-- CorrelationFunctionImpl
    FitterImpl <-- CorrelationFunctionImpl
    FitterImpl <-- LikelihoodImpl
    SourceFunction1D <|-- SourceFunctionImpl
    SourceFunction3D <|-- SourceFunctionImpl
    LikelihoodImpl <|-- ChiSquaredTest
    LikelihoodImpl <|-- LogLikelihoodTest
    class FitterImpl{
        +void Fit()
    }
    class Fitter1D{
        -TH1 data
        +TH1 GetFitResult()
    }
    class Fitter3D{
        -TH3 data
        +TH3 GetFitResult()
    }
    class SourceFunctionImpl{
        +double Value()
    }
    class SourceFunction1D{
        +double Value(radius)
    }
    class SourceFunction3D{
        +double Value(radii1,radii2,radii3)
    }
    class CorrelationFunctionImpl{
    }
    class CorrelationFunction1D{
    }
    class CorrelationFunction3D{
    }
    class InteractionTermImpl{
    }
    class InteractionTermSchrodinger{
    }
    class LikelihoodImpl{
    }
    class ChiSquaredTest{
    }
    class LogLikelihoodTest{
    }

```

## Idea 2 - Template Based

```mermaid
classDiagram

    1D <|-- Dimension
    3D <|-- Dimension
    Fitter~Dimension~ <-- CorrelationFunction~Dimension~
    CorrelationFunction~Dimension~ <-- SourceFunction~Dimension~
    CorrelationFunction~Dimension~ <-- InteractionTerm~Dimension~

    class Dimension{
    }
    class 1D{
    }
    class 3D{
    }
    class Fitter{
        +void Fit()
    }
    class SourceFunction{
        +double Value()
    }
    class CorrelationFunction{
    }
    class InteractionTerm{
    }
```