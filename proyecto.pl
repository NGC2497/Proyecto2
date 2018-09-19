%commands to the interpreter are submitted from stdin input ('show input' box below)
%'halt.' will be automatically appended to stdin input.
%swi-prolog 7.2.3

% Gives list without its last element
removeLast([_], []) :- !.
removeLast([X|Xs], [X|Removed]) :-
    removeLast(Xs, Removed).

% Rename of nth0
getAt(List, I, Res) :-
    nth0(I, List, Res).

% Combines two lists
combine([], L, L) :- !.
combine([H|T], L, [H|M]) :-
    combine(T, L, M).

% Combines two lists and places the new list at last
push(T, L, M) :-
    combine(L, T, M).

% Fills a list with X value N times
fill(_, 0, List) :- not(is_list(List)) -> List = [] ; !.
fill(X, N, List) :-
    % Fills the last list
    M is N - 1,
    fill(X, M, LastList),
    % Last combine
    push([X], LastList, List).

% Concats all the strings in the list
concat([], Str) :- Str = "".
concat(List, Str) :-
    getAt(List, 0, Str0),
    [_|T] = List,
    concat(T, Str1),
    string_concat(Str0, Str1, Str).

% Gets the max of two numbers
max(A, B, Max) :-
A > B -> Max is A ; Max is B.

% Sets the polynomial
polynomial(M_Coeff, M_Deg, Pol) :-
    % Fills list with 0
    fill(0, M_Deg, PreCoeffs),

    % Last value is the coefficient
    push([M_Coeff], PreCoeffs, Coeffs),

    % Gets the degree from the function
    OldPol = [Coeffs, 0],
    p_degree(OldPol, Deg),

    % Returns
    Pol = [Coeffs, Deg].


% Quick access to the coefficients
m_coeffs(Pol, Coeffs) :-
    getAt(Pol, 0, Coeffs).

% Quick access to the degree
m_degree(Pol, Deg) :-
    getAt(Pol, 1, Deg).

% Checks if the polynomial is zero
p_isZero(P, Res) :-
    m_coeffs(P, Coeffs),
    m_degree(P, Deg),
    getAt(Coeffs, 0, Coeff),
    (Deg =:= 0 ->
        (Coeff =:= 0 ->
            Res is 1 ;
            Res is 0) ;
        Res is 0).
    

% Sets the - sign to the polynomial
p_negative([[], _], R) :- R = [[], 0].
p_negative(P, R) :-
    m_coeffs(P, P_Coeffs),
    m_degree(P, P_Deg),
    getAt(P_Coeffs, 0, P_Coeff),
    % Changes sign
    MCoeff is -P_Coeff,
    % Calls for other coefficients
    [_|T] = P_Coeffs,
    p_negative([T, P_Deg], NewR),
    m_coeffs(NewR, NewCoeffs),
    combine([MCoeff], NewCoeffs, R_Coeffs),
    R = [R_Coeffs, P_Deg].


% Gets the degree of a polynomial
p_degree([[], _], _, Deg) :- Deg is 0.
p_degree(Pol, I, Deg) :-
    % Gets coefs
    getAt(Pol, 0, Coeffs),

    % Checks coef
    getAt(Coeffs, 0, Coeff),

    % Checks future coefs
    [_|T] = Coeffs,
    J is I + 1,
    p_degree([T, 0], J, FutureDeg),
    (FutureDeg > I -> Deg is FutureDeg ;
    Coeff \== 0 -> Deg is I ; Deg is 0).
% Just calls the recursive rule
p_degree(Pol, Deg) :-
    p_degree(Pol, 0, Deg).


% Sums two polynomials.
p_sum(P, [[], _], R) :- R = P.
p_sum([[], _], Q, R) :- R = Q.
p_sum(P, Q, R) :-
    m_degree(P, P_Deg),
    m_degree(Q, Q_Deg),
    % Reverts vars if P's degree is higher
    % Gets vars
    m_coeffs(P, P_Coeffs),
    m_coeffs(Q, Q_Coeffs),

    % Gets first coefficients and sums them
    getAt(P_Coeffs, 0, P_Coeff),
    getAt(Q_Coeffs, 0, Q_Coeff),
    R_Coeff is P_Coeff + Q_Coeff,

    % Sums the next coefficients
    [_|PT] = P_Coeffs,
    [_|QT] = Q_Coeffs,
    p_sum([PT, P_Deg], [QT, Q_Deg], NewR),
    getAt(NewR, 0, NewCoeffs),
    % Assigns the new coeffs
    combine([R_Coeff], NewCoeffs, R_Coeffs),
    % Gets new degree
    p_degree([R_Coeffs, 0], R_Deg),
    R = [R_Coeffs, R_Deg].
% Sums all polynomials of the list
p_sum([], R) :- polynomial(0, 0, R).
p_sum(List, R) :-
    % Gets the polynomials
    getAt(List, 0, P),
    [_|T] = List,
    p_sum(T, Q),
    p_sum(P, Q, R).


% Rests P and Q
p_minus(P, Q, R) :-
    p_negative(Q, MQ),
    p_sum(P, MQ, R).


% Multiplies polynomial with scalar
p_scalar([[], _], _, R) :- R = [[], 0].
p_scalar(P, A, R) :-
    %When scalar is 0
    (A =:= 0 ->
        R = [[0], 0] ;
        % Gets vars
        m_coeffs(P, P_Coeffs),

        % Gets first coefficients and sums them
        getAt(P_Coeffs, 0, P_Coeff),
        R_Coeff is P_Coeff * A,

        % Sums the next coefficients
        [_|PT] = P_Coeffs,
        p_scalar([PT, 0], A, NewR),
        getAt(NewR, 0, NewCoeffs),
        % Assigns the new coeffs
        combine([R_Coeff], NewCoeffs, R_Coeffs),
        % Gets new degree
        p_degree([R_Coeffs, 0], R_Deg),
        R = [R_Coeffs, R_Deg]).


% Multiplies two polynomials
p_times([[], _], _, R) :- polynomial(0, 0, R).
p_times(P, Q, R) :-
    % Polynomial 0
    p_isZero(Q, Res),
    (Res =:= 1 ->
        R = [[0], 0] ;
        m_coeffs(P, P_Coeffs),
        % I
        getAt(P_Coeffs, 0, P_Coeff),
        p_scalar(Q, P_Coeff, NewP),
        % J
        [_|T] = P_Coeffs,

        p_times([T, 0], Q, J),
        m_coeffs(J, J_Coeffs),
        combine([0], J_Coeffs, NewQ_Coeffs),
        % Sums
        p_sum(NewP, [NewQ_Coeffs, 0], R)).


% Composes Q in P
p_compose([[], _], _, R) :- polynomial(0, 0, R).
p_compose(P, Q, R) :-
    m_coeffs(P, P_Coeffs),
    getAt(P_Coeffs, 0, P_Coeff),
    % For
    polynomial(P_Coeff, 0, Term),
    [_|T] = P_Coeffs,
    p_compose([T, 0], Q, C),
    p_times(Q, C, NewC),
    p_sum(NewC, Term, R).

% Evaluates polynomial
p_evaluate([[], _], _, Res) :- Res is 0.
p_evaluate(P, X, Res) :-
    m_coeffs(P, Coeffs),
    m_degree(P, Deg),
    getAt(Coeffs, 0, Coeff),
    [_|T] = Coeffs,
    p_evaluate([T, Deg], X, NewRes),
    Res is Coeff + X * NewRes.


% The derivative of a polynomial
p_differentiate(P, D) :-
    m_coeffs(P, P_Coeffs),
    m_degree(P, P_Deg),
    % Degree is 0
    (P_Deg =:= 0 ->
        polynomial(0, 0, D) ;
    % Degree is not 0
    D_Deg is P_Deg - 1,
    getAt(P_Coeffs, P_Deg, P_Coeff),
    D_Coeff is P_Deg * P_Coeff,
    % Calls recursively
    removeLast(P_Coeffs, OldCoeffs),
    p_differentiate([OldCoeffs, D_Deg], NewD),
    m_coeffs(NewD, NewCoeffs),
    (D_Deg \== 0 ->
        combine(NewCoeffs, [D_Coeff], D_Coeffs) ;
        D_Coeffs = [D_Coeff]),
    D = [D_Coeffs, D_Deg]).


% Returns the toString
p_toString(Pol, MaxDeg, Str) :-
    m_coeffs(Pol, Coeffs),
    p_degree(Pol, Deg),
    % Case of degree 0
    (Deg =:= 0 ->
        getAt(Coeffs, 0, Str) ;
    % Higher degree
        getAt(Coeffs, Deg, Coeff),
        % Calls recursive to other coefficients
        removeLast(Coeffs, NewCoeffs),
        % Will create an auxiliary polynomial
        p_degree([NewCoeffs, 0], NewDeg),
        getAt(NewCoeffs, NewDeg, NewCoeff),
        % Ignores if it's zero
        (NewCoeff \== 0 ->
            % Special case when new degree is 0
            (NewDeg =:= 0 ->
                (NewCoeff < 0 ->
                    MNewCoeff is -NewCoeff ;
                    MNewCoeff is NewCoeff),
                p_toString([[MNewCoeff], NewDeg], MaxDeg, NewStr) ;
                p_toString([NewCoeffs, NewDeg], MaxDeg, NewStr)),
            % Sign of the coefficient
            (Coeff < 0 ->
                % Fixes first coefficient
                (Deg =:= MaxDeg ->
                    MCoeff is Coeff ; MCoeff is -Coeff),
                concat([" - ", NewStr], Str0) ;
                (Coeff > 0 ->
                MCoeff is Coeff,
                concat([" + ", NewStr], Str0))) ;
                MCoeff is Coeff,
                Str0 = ""),
        (MCoeff \== 0 ->
            (Deg > 1 ->
            concat([MCoeff, "x^", Deg], Str1) ;
            concat([MCoeff, "x"], Str1)) ; Str1 = ""),
        concat([Str1, Str0], Str)).
p_toString(Pol, Str) :-
    m_degree(Pol, Deg),
    p_toString(Pol, Deg, Str).


% Calls the toString and prints
p_write(Pol) :-
    p_toString(Pol, Str),
    write(Str).

main :-
    polynomial(0, 0, Zero),

    polynomial(4, 3, P1),
    polynomial(3, 2, P2),
    polynomial(1, 0, P3),
    polynomial(2, 1, P4),
    p_sum([P4, P3, P2, P1], P),

    polynomial(3, 2, Q1),
    polynomial(5, 0, Q2),
    p_sum(Q1, Q2, Q),

    p_sum(P, Q, R),
    p_times(P, Q, S),
    p_compose(P, Q, T),
    write("zero(x)     = "), p_write(Zero), nl,
    write("p(x)        = "), p_write(P), nl,
    write("q(x)        = "), p_write(Q), nl,
    write("p(x) + q(x) = "), p_write(R), nl,
    write("p(x) * q(x) = "), p_write(S), nl,
    write("p(q(x))     = "), p_write(T), nl,
    write("0 - p(x)    = "), p_minus(Zero, P, M), p_write(M), nl,
    write("p(3)        = "), p_evaluate(P, 3, E), write(E), nl,
    write("p'(x)       = "), p_differentiate(P, D1), p_write(D1), nl,
    write("p''(x)      = "), p_differentiate(D1, D2), p_write(D2),
    !.
    
:- main.